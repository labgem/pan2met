# Basic profiling of metabolism inference

Install [SnakeViz](https://jiffyclub.github.io/snakeviz/).

Run cProfile
```bash
python3 -m cProfile -o program.prof performance/profile_metabolism_inference.py
```

Launch SnakeViz
```bash
snakeviz program.prof
```
