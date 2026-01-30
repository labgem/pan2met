# Basic profiling of metabolome inference

Install [SnakeViz](https://jiffyclub.github.io/snakeviz/).

Run cProfile
```bash
python3 -m cProfile -o program.prof performance/profile_metabolome_inference.py
```

Launch SnakeViz
```bash
snakeviz program.prof
```
