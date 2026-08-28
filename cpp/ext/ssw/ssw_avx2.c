/*
 *  ssw_avx2.c
 *
 *  AVX2 (256-bit) implementation of striped Smith-Waterman.
 *  Ported from the SSE2 implementation in ssw.c.
 *
 *  This file must be compiled with -mavx2.
 */

#include "ssw_avx2.h"
#include <stdlib.h>
#include <string.h>

#ifdef __GNUC__
#define LIKELY(x) __builtin_expect((x),1)
#define UNLIKELY(x) __builtin_expect((x),0)
#else
#define LIKELY(x) (x)
#define UNLIKELY(x) (x)
#endif

/* Cross-lane byte shift helpers.
 * AVX2 _mm256_slli_si256 shifts within each 128-bit lane independently,
 * so we need explicit cross-lane carry for the striped pattern. */

static inline __m256i avx2_shift_left_1_byte(__m256i v) {
	/* Shift each 128-bit lane left by 1 byte (zeros inserted at byte 0 and byte 16) */
	__m256i shifted = _mm256_slli_si256(v, 1);
	/* Move lower 128-bit lane to upper position, zero the lower lane */
	__m256i carry = _mm256_permute2x128_si256(v, v, 0x08);
	/* Isolate byte[15] from (now in upper lane) */
	carry = _mm256_srli_si256(carry, 15);
	return _mm256_or_si256(shifted, carry);
}

static inline __m256i avx2_shift_left_2_bytes(__m256i v) {
	__m256i shifted = _mm256_slli_si256(v, 2);
	__m256i carry = _mm256_permute2x128_si256(v, v, 0x08);
	carry = _mm256_srli_si256(carry, 14);
	return _mm256_or_si256(shifted, carry);
}

/* Allocate 32-byte aligned memory, zero-initialized */
static inline void* aligned_calloc(size_t count, size_t size) {
	size_t total = count * size;
	/* aligned_alloc requires size to be a multiple of alignment */
	size_t aligned_total = (total + 31) & ~((size_t)31);
	if (aligned_total == 0) aligned_total = 32;
	void* p = aligned_alloc(32, aligned_total);
	if (p) memset(p, 0, aligned_total);
	return p;
}

/* Generate query profile for byte (8-bit) mode.
 * 32 elements per __m256i vector. */
__m256i* qP_avx2_byte(const int8_t* read_num,
                       const int8_t* mat,
                       const int32_t readLen,
                       const int32_t n,
                       uint8_t bias) {

	int32_t segLen = (readLen + 31) / 32;
	__m256i* vProfile = (__m256i*)aligned_calloc(n * segLen, sizeof(__m256i));
	int8_t* t = (int8_t*)vProfile;
	int32_t nt, i, j, segNum;

	for (nt = 0; LIKELY(nt < n); nt++) {
		for (i = 0; i < segLen; i++) {
			j = i;
			for (segNum = 0; LIKELY(segNum < 32); segNum++) {
				*t++ = j >= readLen ? bias : mat[nt * n + read_num[j]] + bias;
				j += segLen;
			}
		}
	}
	return vProfile;
}

/* Generate query profile for word (16-bit) mode.
 * 16 elements per __m256i vector. */
__m256i* qP_avx2_word(const int8_t* read_num,
                       const int8_t* mat,
                       const int32_t readLen,
                       const int32_t n) {

	int32_t segLen = (readLen + 15) / 16;
	__m256i* vProfile = (__m256i*)aligned_calloc(n * segLen, sizeof(__m256i));
	int16_t* t = (int16_t*)vProfile;
	int32_t nt, i, j, segNum;

	for (nt = 0; LIKELY(nt < n); nt++) {
		for (i = 0; i < segLen; i++) {
			j = i;
			for (segNum = 0; LIKELY(segNum < 16); segNum++) {
				*t++ = j >= readLen ? 0 : mat[nt * n + read_num[j]];
				j += segLen;
			}
		}
	}
	return vProfile;
}

/* Striped Smith-Waterman, byte (8-bit) mode, AVX2.
 * Processes 32 query positions in parallel per vector. */
alignment_end_avx2* sw_avx2_byte(const int8_t* ref,
                                  int8_t ref_dir,
                                  int32_t refLen,
                                  int32_t readLen,
                                  const uint8_t weight_gapO,
                                  const uint8_t weight_gapE,
                                  const __m256i* vProfile,
                                  uint8_t terminate,
                                  uint8_t bias,
                                  int32_t maskLen) {

	/* Horizontal max: reduce 32 bytes in vm to scalar m */
	#define max32(m, vm) \
		(vm) = _mm256_max_epu8((vm), _mm256_permute2x128_si256((vm), (vm), 0x01)); \
		(vm) = _mm256_max_epu8((vm), _mm256_srli_si256((vm), 8)); \
		(vm) = _mm256_max_epu8((vm), _mm256_srli_si256((vm), 4)); \
		(vm) = _mm256_max_epu8((vm), _mm256_srli_si256((vm), 2)); \
		(vm) = _mm256_max_epu8((vm), _mm256_srli_si256((vm), 1)); \
		(m) = (uint8_t)_mm256_extract_epi8((vm), 0)

	uint8_t max = 0;
	int32_t end_read = readLen - 1;
	int32_t end_ref = -1;
	int32_t segLen = (readLen + 31) / 32;

	uint8_t* maxColumn = (uint8_t*)calloc(refLen, 1);
	int32_t* end_read_column = (int32_t*)calloc(refLen, sizeof(int32_t));

	__m256i vZero = _mm256_setzero_si256();

	__m256i* pvHStore = (__m256i*)aligned_calloc(segLen, sizeof(__m256i));
	__m256i* pvHLoad  = (__m256i*)aligned_calloc(segLen, sizeof(__m256i));
	__m256i* pvE      = (__m256i*)aligned_calloc(segLen, sizeof(__m256i));
	__m256i* pvHmax   = (__m256i*)aligned_calloc(segLen, sizeof(__m256i));

	int32_t i, j, k;
	__m256i vGapO = _mm256_set1_epi8(weight_gapO);
	__m256i vGapE = _mm256_set1_epi8(weight_gapE);
	__m256i vBias = _mm256_set1_epi8(bias);

	__m256i vMaxScore = vZero;
	__m256i vMaxMark = vZero;
	__m256i vTemp;
	int32_t edge, begin = 0, end = refLen, step = 1;

	if (ref_dir == 1) {
		begin = refLen - 1;
		end = -1;
		step = -1;
	}
	for (i = begin; LIKELY(i != end); i += step) {
		int32_t cmp;
		__m256i e, vF = vZero, vMaxColumn = vZero;

		__m256i vH = pvHStore[segLen - 1];
		vH = avx2_shift_left_1_byte(vH);
		const __m256i* vP = vProfile + ref[i] * segLen;

		__m256i* pv = pvHLoad;
		pvHLoad = pvHStore;
		pvHStore = pv;

		/* Inner loop: process query segments */
		for (j = 0; LIKELY(j < segLen); ++j) {
			vH = _mm256_adds_epu8(vH, _mm256_load_si256(vP + j));
			vH = _mm256_subs_epu8(vH, vBias);

			e = _mm256_load_si256(pvE + j);
			vH = _mm256_max_epu8(vH, e);
			vH = _mm256_max_epu8(vH, vF);
			vMaxColumn = _mm256_max_epu8(vMaxColumn, vH);

			_mm256_store_si256(pvHStore + j, vH);

			vH = _mm256_subs_epu8(vH, vGapO);
			e = _mm256_subs_epu8(e, vGapE);
			e = _mm256_max_epu8(e, vH);
			_mm256_store_si256(pvE + j, e);

			vF = _mm256_subs_epu8(vF, vGapE);
			vF = _mm256_max_epu8(vF, vH);

			vH = _mm256_load_si256(pvHLoad + j);
		}

		/* Lazy_F loop */
		for (k = 0; LIKELY(k < 32); ++k) {
			vF = avx2_shift_left_1_byte(vF);
			for (j = 0; LIKELY(j < segLen); ++j) {
				vH = _mm256_load_si256(pvHStore + j);
				vH = _mm256_max_epu8(vH, vF);
				vMaxColumn = _mm256_max_epu8(vMaxColumn, vH);
				_mm256_store_si256(pvHStore + j, vH);
				vH = _mm256_subs_epu8(vH, vGapO);
				vF = _mm256_subs_epu8(vF, vGapE);
				vTemp = _mm256_subs_epu8(vF, vH);
				vTemp = _mm256_cmpeq_epi8(vTemp, vZero);
				if (UNLIKELY(_mm256_movemask_epi8(vTemp) == (int32_t)0xffffffff)) goto end;
			}
		}

end:
		vMaxScore = _mm256_max_epu8(vMaxScore, vMaxColumn);
		vTemp = _mm256_cmpeq_epi8(vMaxMark, vMaxScore);
		cmp = _mm256_movemask_epi8(vTemp);
		if (cmp != (int32_t)0xffffffff) {
			uint8_t temp;
			vMaxMark = vMaxScore;
			max32(temp, vMaxScore);
			vMaxScore = vMaxMark;

			if (LIKELY(temp > max)) {
				max = temp;
				if (max + bias >= 255) break;
				end_ref = i;
				for (j = 0; LIKELY(j < segLen); ++j) pvHmax[j] = pvHStore[j];
			}
		}

		/* Record the max score of current column. */
		max32(maxColumn[i], vMaxColumn);
		if (maxColumn[i] == terminate) break;
	}

	/* Trace the alignment ending position on read. */
	uint8_t *t = (uint8_t*)pvHmax;
	int32_t column_len = segLen * 32;
	for (i = 0; LIKELY(i < column_len); ++i, ++t) {
		int32_t temp;
		if (*t == max) {
			temp = i / 32 + i % 32 * segLen;
			if (temp < end_read) end_read = temp;
		}
	}

	free(pvHmax);
	free(pvE);
	free(pvHLoad);
	free(pvHStore);

	alignment_end_avx2* bests = (alignment_end_avx2*)calloc(2, sizeof(alignment_end_avx2));
	bests[0].score = max + bias >= 255 ? 255 : max;
	bests[0].ref = end_ref;
	bests[0].read = end_read;

	bests[1].score = 0;
	bests[1].ref = 0;
	bests[1].read = 0;

	edge = (end_ref - maskLen) > 0 ? (end_ref - maskLen) : 0;
	for (i = 0; i < edge; i++) {
		if (maxColumn[i] > bests[1].score) {
			bests[1].score = maxColumn[i];
			bests[1].ref = i;
		}
	}
	edge = (end_ref + maskLen) > refLen ? refLen : (end_ref + maskLen);
	for (i = edge + 1; i < refLen; i++) {
		if (maxColumn[i] > bests[1].score) {
			bests[1].score = maxColumn[i];
			bests[1].ref = i;
		}
	}

	free(maxColumn);
	free(end_read_column);
	return bests;

	#undef max32
}

/* Striped Smith-Waterman, word (16-bit) mode, AVX2.
 * Processes 16 query positions in parallel per vector. */
alignment_end_avx2* sw_avx2_word(const int8_t* ref,
                                  int8_t ref_dir,
                                  int32_t refLen,
                                  int32_t readLen,
                                  const uint8_t weight_gapO,
                                  const uint8_t weight_gapE,
                                  const __m256i* vProfile,
                                  uint16_t terminate,
                                  int32_t maskLen) {

	/* Horizontal max: reduce 16 words in vm to scalar m */
	#define max16_avx2(m, vm) \
		(vm) = _mm256_max_epi16((vm), _mm256_permute2x128_si256((vm), (vm), 0x01)); \
		(vm) = _mm256_max_epi16((vm), _mm256_srli_si256((vm), 8)); \
		(vm) = _mm256_max_epi16((vm), _mm256_srli_si256((vm), 4)); \
		(vm) = _mm256_max_epi16((vm), _mm256_srli_si256((vm), 2)); \
		(m) = (uint16_t)_mm256_extract_epi16((vm), 0)

	uint16_t max = 0;
	int32_t end_read = readLen - 1;
	int32_t end_ref = 0;
	int32_t segLen = (readLen + 15) / 16;

	uint16_t* maxColumn = (uint16_t*)calloc(refLen, sizeof(uint16_t));
	int32_t* end_read_column = (int32_t*)calloc(refLen, sizeof(int32_t));

	__m256i vZero = _mm256_setzero_si256();

	__m256i* pvHStore = (__m256i*)aligned_calloc(segLen, sizeof(__m256i));
	__m256i* pvHLoad  = (__m256i*)aligned_calloc(segLen, sizeof(__m256i));
	__m256i* pvE      = (__m256i*)aligned_calloc(segLen, sizeof(__m256i));
	__m256i* pvHmax   = (__m256i*)aligned_calloc(segLen, sizeof(__m256i));

	int32_t i, j, k;
	__m256i vGapO = _mm256_set1_epi16(weight_gapO);
	__m256i vGapE = _mm256_set1_epi16(weight_gapE);

	__m256i vMaxScore = vZero;
	__m256i vMaxMark = vZero;
	__m256i vTemp;
	int32_t edge, begin = 0, end = refLen, step = 1;

	if (ref_dir == 1) {
		begin = refLen - 1;
		end = -1;
		step = -1;
	}
	for (i = begin; LIKELY(i != end); i += step) {
		int32_t cmp;
		__m256i e, vF = vZero;

		__m256i vH = pvHStore[segLen - 1];
		vH = avx2_shift_left_2_bytes(vH);

		__m256i* pv = pvHLoad;

		__m256i vMaxColumn = vZero;

		const __m256i* vP = vProfile + ref[i] * segLen;
		pvHLoad = pvHStore;
		pvHStore = pv;

		/* Inner loop: process query segments */
		for (j = 0; LIKELY(j < segLen); j++) {
			vH = _mm256_adds_epi16(vH, _mm256_load_si256(vP + j));

			e = _mm256_load_si256(pvE + j);
			vH = _mm256_max_epi16(vH, e);
			vH = _mm256_max_epi16(vH, vF);
			vMaxColumn = _mm256_max_epi16(vMaxColumn, vH);

			_mm256_store_si256(pvHStore + j, vH);

			vH = _mm256_subs_epu16(vH, vGapO);
			e = _mm256_subs_epu16(e, vGapE);
			e = _mm256_max_epi16(e, vH);
			_mm256_store_si256(pvE + j, e);

			vF = _mm256_subs_epu16(vF, vGapE);
			vF = _mm256_max_epi16(vF, vH);

			vH = _mm256_load_si256(pvHLoad + j);
		}

		/* Lazy_F loop */
		for (k = 0; LIKELY(k < 16); ++k) {
			vF = avx2_shift_left_2_bytes(vF);
			for (j = 0; LIKELY(j < segLen); ++j) {
				vH = _mm256_load_si256(pvHStore + j);
				vH = _mm256_max_epi16(vH, vF);
				vMaxColumn = _mm256_max_epi16(vMaxColumn, vH);
				_mm256_store_si256(pvHStore + j, vH);
				vH = _mm256_subs_epu16(vH, vGapO);
				vF = _mm256_subs_epu16(vF, vGapE);
				if (UNLIKELY(!_mm256_movemask_epi8(_mm256_cmpgt_epi16(vF, vH)))) goto end;
			}
		}

end:
		vMaxScore = _mm256_max_epi16(vMaxScore, vMaxColumn);
		vTemp = _mm256_cmpeq_epi16(vMaxMark, vMaxScore);
		cmp = _mm256_movemask_epi8(vTemp);
		if (cmp != (int32_t)0xffffffff) {
			uint16_t temp;
			vMaxMark = vMaxScore;
			max16_avx2(temp, vMaxScore);
			vMaxScore = vMaxMark;

			if (LIKELY(temp > max)) {
				max = temp;
				end_ref = i;
				for (j = 0; LIKELY(j < segLen); ++j) pvHmax[j] = pvHStore[j];
			}
		}

		/* Record the max score of current column. */
		max16_avx2(maxColumn[i], vMaxColumn);
		if (maxColumn[i] == terminate) break;
	}

	/* Trace the alignment ending position on read. */
	uint16_t *t = (uint16_t*)pvHmax;
	int32_t column_len = segLen * 16;
	for (i = 0; LIKELY(i < column_len); ++i, ++t) {
		int32_t temp;
		if (*t == max) {
			temp = i / 16 + i % 16 * segLen;
			if (temp < end_read) end_read = temp;
		}
	}

	free(pvHmax);
	free(pvE);
	free(pvHLoad);
	free(pvHStore);

	alignment_end_avx2* bests = (alignment_end_avx2*)calloc(2, sizeof(alignment_end_avx2));
	bests[0].score = max;
	bests[0].ref = end_ref;
	bests[0].read = end_read;

	bests[1].score = 0;
	bests[1].ref = 0;
	bests[1].read = 0;

	edge = (end_ref - maskLen) > 0 ? (end_ref - maskLen) : 0;
	for (i = 0; i < edge; i++) {
		if (maxColumn[i] > bests[1].score) {
			bests[1].score = maxColumn[i];
			bests[1].ref = i;
		}
	}
	edge = (end_ref + maskLen) > refLen ? refLen : (end_ref + maskLen);
	for (i = edge; i < refLen; i++) {
		if (maxColumn[i] > bests[1].score) {
			bests[1].score = maxColumn[i];
			bests[1].ref = i;
		}
	}

	free(maxColumn);
	free(end_read_column);
	return bests;

	#undef max16_avx2
}
