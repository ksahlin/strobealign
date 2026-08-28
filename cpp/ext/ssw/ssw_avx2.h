#ifndef SSW_AVX2_H
#define SSW_AVX2_H

#include <stdint.h>
#include <immintrin.h>

typedef struct {
	uint16_t score;
	int32_t ref;
	int32_t read;
} alignment_end_avx2;

__m256i* qP_avx2_byte(const int8_t* read_num, const int8_t* mat,
                       int32_t readLen, int32_t n, uint8_t bias);

__m256i* qP_avx2_word(const int8_t* read_num, const int8_t* mat,
                       int32_t readLen, int32_t n);

alignment_end_avx2* sw_avx2_byte(const int8_t* ref, int8_t ref_dir,
                                  int32_t refLen, int32_t readLen,
                                  uint8_t weight_gapO, uint8_t weight_gapE,
                                  const __m256i* vProfile, uint8_t terminate,
                                  uint8_t bias, int32_t maskLen);

alignment_end_avx2* sw_avx2_word(const int8_t* ref, int8_t ref_dir,
                                  int32_t refLen, int32_t readLen,
                                  uint8_t weight_gapO, uint8_t weight_gapE,
                                  const __m256i* vProfile, uint16_t terminate,
                                  int32_t maskLen);

#endif
