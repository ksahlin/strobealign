#include "insertsizedistribution.hpp"
#include <cmath>
#include <iostream>


/* Add a new observation */
void InsertSizeDistribution::update(int dist) {
    if (dist >= 2000) {
        return;
    }
    const float e = dist - mu;
    mu += e / sample_size; // (1.0/(sample_size +1.0)) * (sample_size*mu + d);
    SSE += e * (dist - mu);
    if (sample_size > 1) {
        //d < 1000 ? ((sample_size +1.0)/sample_size) * ( (V*sample_size/(sample_size +1)) + ((mu-d)*(mu-d))/sample_size ) : V;
        V = SSE / (sample_size - 1.0);
    } else {
        V = SSE;
    }
    sigma = std::sqrt(V);
    sample_size = sample_size + 1.0;
    if (mu < 0) {
        std::cerr << "mu negative, mu: " << mu << " sigma: " << sigma << " SSE: " << SSE << " sample size: " << sample_size << std::endl;
    }
    if (SSE < 0) {
        std::cerr << "SSE negative, mu: " << mu << " sigma: " << sigma << " SSE: " << SSE << " sample size: " << sample_size << std::endl;
    }
}
