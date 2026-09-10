# Contributing

Contributions are welcome, and they are greatly appreciated!
Every little bit helps, and credit will always be given.

## Environment setup

Fork and clone the repository:

```bash
git clone https://github.com/labgem/pan2met.git
cd pan2met
```

We use [pixi](prefix.dev) package management tool to manage our dependencies and virtual environment.

You can run `pixi shell` from the repository root and you should have the dependencies installed.

You can run the tests with `python3 -m pytest`.

## Development

As usual:

1. create a new branch: `git checkout -b feature-or-bugfix-name`
2. edit the code and/or the documentation

If you updated the documentation or the project dependencies:

1. run `cd docs/ && make html && python3 -m http.server 8000`
2. go to http://localhost:8000 and check that everything looks good
