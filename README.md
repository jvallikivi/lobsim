# Simple Experimental Limit Order Book Simulation

### Requirements:
   * Rust
   * Python (plotly, cffi)
   
### Run:
For a simple demo run:

```
cargo build --release
python lobsim.py
```

Modify values inside the [file](lobsim.py) to experiment with different values and to
save results.

### Notes:

The distributions for this simulation have been taken from LOBSTER samples. 
See [https://lobsterdata.com/](https://lobsterdata.com/) for more.

This model is based on the queue-reactive model for simulating order book data. 
See more:

Huang, W., Lehalle, C.A. and Rosenbaum, M., 2015. Simulating and analyzing order book data: The queue-reactive model. Journal of the American Statistical Association, 110(509), pp.107-122.

[Link](https://arxiv.org/abs/1312.0563)
