// Compile and run with:

#include <cudaq.h>

struct kernel {
  void operator()(const int N, const float time_param, const int steps) __qpu__ {
    cudaq::qreg q(N);

    int k0 = 0, k1 = 0;

    for (int i = 0; i < steps; ++i) {
      
      // Prepare the initial state by superposition
      h(q[1]);

      // Problem Hamiltonian: ZY
      rx(M_PI/2, q[1]);
      x<cudaq::ctrl>(q[1], q[0]);
      rz(2*time_param, q[0]);
      x<cudaq::ctrl>(q[1], q[0]);
      rx(-M_PI/2, q[1]);


      // Mid circuit measurement
      auto c = mz(q[1]);

      // Updating k0 and k1
      int _k0 = k0 + (int)(0 == c);
      int _k1 = k1 + (int)(1 == c);
      float val = (1.0*(_k1-_k0))/(_k0+_k1+0.0001);
      float y = val*(1.0+val*val*(1.0/6 + val*val*(3.0/(2*4*5) + val*val*((1.0*3*5)/(2*4*6*7))))); // arcsin 
      float x_max = 0.5*y;

      if(x_max < 0){
        k0 += (int)(0 == c);
        k1 += (int)(1 == c);
      }
      else{
        k0 = 0;
        k1 = 0;
      }

      // Apply conditional unitary
      val = (1.0*(k1-k0))/(k0+k1+0.0001);
      y = val*(1.0+val*val*(1.0/6 + val*val*(3.0/(2*4*5) + val*val*((1.0*3*5)/(2*4*6*7))))); // arcsin 
      x_max = 0.5*y;

      if(x_max >= 0){
        x(q[0]);
      }

      // arcsin implementation. ref: https://dsp.stackexchange.com/questions/25770/looking-for-an-arcsin-algorithm
      // unable to use math.h asin. Error: undefined reference to `asin'

      reset(q[1]);   
    }
  }
};

int main() {

  using namespace cudaq::spin;

  // Problem parameters
  const int n_qubits = 2;
  const int nShots = 100;
  const float epsilon = 0.2;
  const int steps = 100;


  // Sample
  auto counts = cudaq::sample(nShots, kernel{}, n_qubits, epsilon, steps);

  // Dump the states and the counts
  // counts.dump();

  // Fine grain access to the bits and counts
  for (auto &[bits, count] : counts) {
    printf("Observed: %s, %lu\n", bits.data(), count);
  }

}

