/// Estimator for a normal distribution, used for insert sizes.
#[derive(Debug)]
pub struct InsertSizeDistribution {
    pub sample_size: usize,
    pub mu: f32,
    pub sigma: f32,
    V: f32,
    SSE: f32,
}

impl Default for InsertSizeDistribution {
    fn default() -> Self {
        InsertSizeDistribution {
            sample_size: 0,
            mu: 300.0,
            sigma: 100.0,
            V: 10_000.0,
            SSE: 10_000.0,
        }
    }
}

impl InsertSizeDistribution {
    pub fn new() -> Self {
        InsertSizeDistribution::default()
    }

    /// Add a new observation
    pub fn update(&mut self, insert_size: usize) {
        if insert_size >= 2000 {
            return;
        }
        if self.sample_size > 0 {
            let insert_size = insert_size as f32;
            let e = insert_size - self.mu;
            self.mu += e / self.sample_size as f32;
            self.SSE += e * (insert_size - self.mu);
            if self.sample_size > 1 {
                self.V = self.SSE / (self.sample_size as f32 - 1.0);
            } else {
                self.V = self.SSE;
            }
            self.sigma = self.V.sqrt();
        }
        self.sample_size += 1;
        /*
        TODO
        if self.mu < 0 {
            std::cerr << "mu negative, mu: " << mu << " sigma: " << sigma << " SSE: " << SSE << " sample size: " << sample_size << std::endl;
        }
        if self.SSE < 0 {
            std::cerr << "SSE negative, mu: " << mu << " sigma: " << sigma << " SSE: " << SSE << " sample size: " << sample_size << std::endl;
        }
        */
    }
}
