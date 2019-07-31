# Running Benchmarks

We use `criterion` crate for performing benchmarks.
Benchmarks can be run via `cargo bench`. To run a specific
set of benchmarks:

```
cargo bench --bench delaunay -- <filter>
```

For instance `f64` to run only floating-point tests, or
`i64` for only integer tests.

# Delaunay Benchmarks

We run the benchmarks on all the available kernel,
locate-strategy combinations, and at varying sizes.

## Insertion

We benchmark the time taken to insert n random values, where
n quadruples from 1 to 16384 (8 data points in total).

### Float Kernels

```
insert/uniform_f64/float_kernel/tree_locate/1           1.00   247.7±34.48ns        ? B/sec
insert/uniform_f64/float_kernel/tree_locate/4           1.00  1359.9±43.54ns        ? B/sec
insert/uniform_f64/float_kernel/tree_locate/16          1.00     10.4±0.25µs        ? B/sec
insert/uniform_f64/float_kernel/tree_locate/64          1.00     73.5±0.84µs        ? B/sec
insert/uniform_f64/float_kernel/tree_locate/256         1.00    402.8±1.16µs        ? B/sec
insert/uniform_f64/float_kernel/tree_locate/1024        1.00      2.0±0.01ms        ? B/sec
insert/uniform_f64/float_kernel/tree_locate/4096        1.00      9.3±0.07ms        ? B/sec
insert/uniform_f64/float_kernel/tree_locate/16384       1.00     43.3±0.56ms        ? B/sec
insert/uniform_f64/float_kernel/walk_locate/1           1.00    116.7±1.10ns        ? B/sec
insert/uniform_f64/float_kernel/walk_locate/4           1.00   982.8±14.96ns        ? B/sec
insert/uniform_f64/float_kernel/walk_locate/16          1.00      4.9±0.26µs        ? B/sec
insert/uniform_f64/float_kernel/walk_locate/64          1.00     22.8±0.61µs        ? B/sec
insert/uniform_f64/float_kernel/walk_locate/256         1.00    171.0±0.62µs        ? B/sec
insert/uniform_f64/float_kernel/walk_locate/1024        1.00   1084.2±1.40µs        ? B/sec
insert/uniform_f64/float_kernel/walk_locate/4096        1.00     10.3±0.03ms        ? B/sec
insert/uniform_f64/float_kernel/walk_locate/16384       1.00     91.8±0.44ms        ? B/sec
insert/uniform_f64/trivial_kernel/tree_locate/1         1.00   290.5±26.25ns        ? B/sec
insert/uniform_f64/trivial_kernel/tree_locate/4         1.00  1348.2±50.04ns        ? B/sec
insert/uniform_f64/trivial_kernel/tree_locate/16        1.00     10.1±0.22µs        ? B/sec
insert/uniform_f64/trivial_kernel/tree_locate/64        1.00     68.3±1.07µs        ? B/sec
insert/uniform_f64/trivial_kernel/tree_locate/256       1.00    374.8±1.14µs        ? B/sec
insert/uniform_f64/trivial_kernel/tree_locate/1024      1.00   1891.8±7.48µs        ? B/sec
insert/uniform_f64/trivial_kernel/tree_locate/4096      1.00      8.8±0.08ms        ? B/sec
insert/uniform_f64/trivial_kernel/tree_locate/16384     1.00     41.0±0.58ms        ? B/sec
insert/uniform_f64/trivial_kernel/walk_locate/1         1.00    115.5±0.84ns        ? B/sec
insert/uniform_f64/trivial_kernel/walk_locate/4         1.00   967.5±17.72ns        ? B/sec
insert/uniform_f64/trivial_kernel/walk_locate/16        1.00      4.5±0.30µs        ? B/sec
insert/uniform_f64/trivial_kernel/walk_locate/64        1.00     19.1±0.80µs        ? B/sec
insert/uniform_f64/trivial_kernel/walk_locate/256       1.00    123.0±0.73µs        ? B/sec
insert/uniform_f64/trivial_kernel/walk_locate/1024      1.00    802.2±3.29µs        ? B/sec
insert/uniform_f64/trivial_kernel/walk_locate/4096      1.00      7.8±0.04ms        ? B/sec
insert/uniform_f64/trivial_kernel/walk_locate/16384     1.00     73.3±0.58ms        ? B/sec
```


### Integer Kernels

```
insert/uniform_i64/adaptive_kernel/tree_locate/1        1.00   264.8±35.94ns        ? B/sec
insert/uniform_i64/adaptive_kernel/tree_locate/4        1.00      2.5±0.37µs        ? B/sec
insert/uniform_i64/adaptive_kernel/tree_locate/16       1.00     63.4±5.57µs        ? B/sec
insert/uniform_i64/adaptive_kernel/tree_locate/64       1.00   470.6±22.35µs        ? B/sec
insert/uniform_i64/adaptive_kernel/tree_locate/256      1.00      2.4±0.05ms        ? B/sec
insert/uniform_i64/adaptive_kernel/tree_locate/1024     1.00     10.8±0.11ms        ? B/sec
insert/uniform_i64/adaptive_kernel/tree_locate/4096     1.00     46.2±0.34ms        ? B/sec
insert/uniform_i64/adaptive_kernel/tree_locate/16384    1.00    196.4±1.26ms        ? B/sec
insert/uniform_i64/adaptive_kernel/walk_locate/1        1.00    114.2±1.02ns        ? B/sec
insert/uniform_i64/adaptive_kernel/walk_locate/4        1.00      2.2±0.46µs        ? B/sec
insert/uniform_i64/adaptive_kernel/walk_locate/16       1.00     54.7±4.88µs        ? B/sec
insert/uniform_i64/adaptive_kernel/walk_locate/64       1.00   394.6±14.43µs        ? B/sec
insert/uniform_i64/adaptive_kernel/walk_locate/256      1.00      2.0±0.03ms        ? B/sec
insert/uniform_i64/adaptive_kernel/walk_locate/1024     1.00      9.2±0.08ms        ? B/sec
insert/uniform_i64/adaptive_kernel/walk_locate/4096     1.00     42.9±0.25ms        ? B/sec
insert/uniform_i64/adaptive_kernel/walk_locate/16384    1.00    217.8±0.93ms        ? B/sec
insert/uniform_i64/trivial_kernel/tree_locate/1         1.00   305.3±32.50ns        ? B/sec
insert/uniform_i64/trivial_kernel/tree_locate/4         1.00  1381.6±57.60ns        ? B/sec
insert/uniform_i64/trivial_kernel/walk_locate/1         1.00     81.0±1.87ns        ? B/sec
insert/uniform_i64/trivial_kernel/walk_locate/4         1.00   743.0±23.15ns        ? B/sec
```

## Locate Queries

We benchmark the time taken to locate a random point, given
a tree with n random points. Again, n quadruples from 1 to
16384 (8 data points in total).

### Float Kernels

```
locate/uniform_f64/float_kernel/tree_locate/1           1.00    513.1±0.67ns        ? B/sec
locate/uniform_f64/float_kernel/tree_locate/4           1.00    585.4±0.49ns        ? B/sec
locate/uniform_f64/float_kernel/tree_locate/16          1.00    642.2±0.61ns        ? B/sec
locate/uniform_f64/float_kernel/tree_locate/64          1.00    690.5±0.71ns        ? B/sec
locate/uniform_f64/float_kernel/tree_locate/256         1.00    763.1±0.60ns        ? B/sec
locate/uniform_f64/float_kernel/tree_locate/1024        1.00    850.2±2.35ns        ? B/sec
locate/uniform_f64/float_kernel/tree_locate/4096        1.00    975.9±3.17ns        ? B/sec
locate/uniform_f64/float_kernel/tree_locate/16384       1.00   1154.4±4.17ns        ? B/sec
locate/uniform_f64/float_kernel/walk_locate/1           1.00    514.2±0.27ns        ? B/sec
locate/uniform_f64/float_kernel/walk_locate/4           1.00    569.6±0.57ns        ? B/sec
locate/uniform_f64/float_kernel/walk_locate/16          1.00    629.0±0.39ns        ? B/sec
locate/uniform_f64/float_kernel/walk_locate/64          1.00    811.2±0.66ns        ? B/sec
locate/uniform_f64/float_kernel/walk_locate/256         1.00    841.4±0.56ns        ? B/sec
locate/uniform_f64/float_kernel/walk_locate/1024        1.00   1010.0±1.12ns        ? B/sec
locate/uniform_f64/float_kernel/walk_locate/4096        1.00      3.0±0.00µs        ? B/sec
locate/uniform_f64/float_kernel/walk_locate/16384       1.00      3.8±0.01µs        ? B/sec
locate/uniform_f64/trivial_kernel/tree_locate/1         1.00    512.6±0.40ns        ? B/sec
locate/uniform_f64/trivial_kernel/tree_locate/4         1.00    580.9±0.49ns        ? B/sec
locate/uniform_f64/trivial_kernel/tree_locate/16        1.00    637.1±2.80ns        ? B/sec
locate/uniform_f64/trivial_kernel/tree_locate/64        1.00    681.5±1.24ns        ? B/sec
locate/uniform_f64/trivial_kernel/tree_locate/256       1.00    736.6±1.32ns        ? B/sec
locate/uniform_f64/trivial_kernel/tree_locate/1024      1.00    816.9±2.19ns        ? B/sec
locate/uniform_f64/trivial_kernel/tree_locate/4096      1.00    937.1±6.65ns        ? B/sec
locate/uniform_f64/trivial_kernel/tree_locate/16384     1.00   1136.4±1.70ns        ? B/sec
locate/uniform_f64/trivial_kernel/walk_locate/1         1.00    509.4±0.43ns        ? B/sec
locate/uniform_f64/trivial_kernel/walk_locate/4         1.00    561.4±0.90ns        ? B/sec
locate/uniform_f64/trivial_kernel/walk_locate/16        1.00    592.4±0.33ns        ? B/sec
locate/uniform_f64/trivial_kernel/walk_locate/64        1.00    704.7±0.22ns        ? B/sec
locate/uniform_f64/trivial_kernel/walk_locate/256       1.00    723.9±0.31ns        ? B/sec
locate/uniform_f64/trivial_kernel/walk_locate/1024      1.00    843.0±0.74ns        ? B/sec
locate/uniform_f64/trivial_kernel/walk_locate/4096      1.00      2.3±0.00µs        ? B/sec
locate/uniform_f64/trivial_kernel/walk_locate/16384     1.00      3.1±0.01µs        ? B/sec
```

### Integer Kernels

```
locate/uniform_i64/adaptive_kernel/tree_locate/1        1.00    509.7±0.36ns        ? B/sec
locate/uniform_i64/adaptive_kernel/tree_locate/4        1.00    588.9±0.47ns        ? B/sec
locate/uniform_i64/adaptive_kernel/walk_locate/1        1.00    505.0±1.25ns        ? B/sec
locate/uniform_i64/adaptive_kernel/walk_locate/4        1.00    548.2±0.31ns        ? B/sec
locate/uniform_i64/trivial_kernel/walk_locate/1         1.00    510.2±0.30ns        ? B/sec
```
