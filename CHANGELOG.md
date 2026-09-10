# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- a metabolic pathway prediction algorithm
    - a pathway score to estimate a confidence in the pathway prediction
- a method to identify predicted chimeric pathway
- using the genomic context to enhance the annotation of enzyme catalysis using a transitive closure on gene node in the pangenome graph
- a method to identify a set of minimal monomer that could explain a set of observed catalyzed reactions
