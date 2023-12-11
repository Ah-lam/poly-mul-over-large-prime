# poly-mul-over-large-prime
This sample code package is an implementation of polynomial multiplication under large prime numbers (Q=524190235384903211525979137). 

The multiplication of these polynomials utilizes the Number Theoretic Transform (NTT) algorithm, capable of handling up to 1024 polynomials represented in the form of (a + bx) during the multiplication process.

## License

This project is licensed under the Apache-2.0 License.


BUILD
-----

```
make main_test
```

To run

`./main_test`



Performance measurements
------------------------
The performance measurements are reported in processor cycles (per single core). The results are obtained using the following methodology. Measured function runs one time (warm-up), followed by 100000 iterations that were clocked and averaged. 

To run the benchmarking

`./bench`