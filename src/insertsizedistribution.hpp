#ifndef STROBEALIGN_INSERTSIZEDISTRIBUTION_HPP
#define STROBEALIGN_INSERTSIZEDISTRIBUTION_HPP

/* Estimator for a normal distribution, used for insert sizes.
 */
class InsertSizeDistribution {
public:
    float sample_size = 1;
    float mu = 300;
    float sigma = 100;
    float V = 10000;
    float SSE = 10000;

    // Add a new observation
    void update(int dist);
};

#endif
