# Reinforcement Learning Under Uncertainty: A Replication Study

This repository contains a direct computational replication of the study by Ez-zizi et al. (2023), which investigated human learning under conditions of perceptual (state) and reward uncertainty. The original work evaluated how well different reinforcement learning models explain human behavior across these conditions.

## Overview

In their study, Ez-zizi et al. (2023) examined how people adapt to:
- **State uncertainty** (ambiguous sensory input)
- **Reward uncertainty** (noisy outcome feedback)

Participants performed reversal learning tasks, and three reinforcement learning models were fit to the data:
1. Simple RL (WTA-RL)
2. Bayesian RL (PW-RL)
3. Sampling-based models

Our goal was to reproduce their behavioral and computational findings using the original data and code.

## Objectives

- Independently validate the findings of Ez-zizi et al. (2023)
- Evaluate model robustness under limited computational resources
- Promote reproducibility and transparency in computational cognitive science

## Replicated Models

- Simple Reinforcement Learning (WTA-RL)
- Bayesian Reinforcement Learning (PW-RL)
- Sampling-based Decision Models

Model performance was compared using AIC scores, and learning dynamics were visualized across conditions.

## Methods

The original dataset and analysis scripts were obtained from the open-access repository referenced in the original paper. We executed all simulations and model fitting using MATLAB and R. Due to limitations in computing power, we reduced the number of agents and iterations in some simulations.

## Repository

Clone this repository using:

```
gh repo clone Yuksel-I-Ozturk-M-2025/replication-rl-uncertainty
```

The project directory includes:
- `/Code`: Original model scripts and preprocessing
- `/Results`: Reproduction results and figures
- `/Docs`: Final paper and supplementary material

## Project Report

The full paper is available on Overleaf:  
[https://www.overleaf.com/read/szqcxwvbrcmb](https://www.overleaf.com/read/szqcxwvbrcmb)

## Students

- İbrahim Yüksel  
- Miraç Öztürk  
Boğaziçi University, Department of Computer Engineering  
CMPE 489 – Cognitive Science, Spring 2025
