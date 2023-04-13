# LEARN
Likelihood Event Analysis of Reactor Neutrinos

### Dependencies
ROOT must be installed and sourced before compiling or using. Developed using ROOT 6.18

### Compile
Compile:
```bash
make
```

Clean:
```bash
make clean
```

Source and export path:
```bash
source LEARN.sh
```

### Use

To see input options:
```bash
learn
```
To print analysis steps:
```bash
learn --help
```

There are several stages to this analysis

1) Creating PDFs:
```bash
learn --pdf -f [input file] -c [component]
```

2) Scaling and merging PDFs:
```bash
learn --merge -s [signal]
```

3) Creating likelihoods:
```bash
learn --like -f [input file] -c [component]
```

4) Evaluating likelihoods and rates:
```bash
learn --eval -f [input file] -c [component]
```

5) Data reduction and fast neutron handling:
Not currently integrated into LEARN

6) Evaluating dwell times (energy cut):
```bash
learn --energy -s [signal]
```

7) Evaluating dwell times (no energy cut):
```bash
learn --veto -s [signal]
```

Event rates and components are defined in 'rates.csv', which needs to be in the current working directory. Example files for different detector configurations are found in this repository.
