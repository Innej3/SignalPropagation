# Signal-Propagation
## Summary

The paper's temporal distance predicts signal arrival times well for regulatory dynamics, but fails for SIS epidemics on scale-free networks, and fails completely when used to locate the source of a perturbation.

I replicated the core results of Hens, Harush, Haber, Cohen and Barzel, *Spatiotemporal signal propagation in complex networks* ([Nature Physics 15, 403–412, 2019](https://arxiv.org/abs/1801.08854)), for three dynamics (two gene-regulation models and SIS) on Erdős–Rényi and scale-free networks. The local scaling of response time with degree reproduces for all three. The collapse of arrival times onto the temporal distance holds for the regulatory models, including the hub-bottlenecked case where hop count is nearly useless (correlation 0.91 for temporal distance vs 0.25 for hops).

For SIS on scale-free networks the picture flips: hop count predicts arrival times better than temporal distance (0.80 vs 0.63). Running the problem backwards makes the gap stark. With 20% of nodes observed, a hop-count estimator finds the true source 89% of the time, the temporal-distance estimator about 1%, and an oracle using the true dynamics 100%. The arrival times carry the information; the paper's metric does not extract it.

Code: two Jupyter notebooks, a replication notebook covering ℝ₁, ℝ₂ and SIS, and an SIS and retrodiction notebook.

## Background

The paper asks: if you perturb node j, when does node i feel it? Network topology alone cannot answer this, because the same graph running different dynamics spreads signals in very different ways.

All dynamics studied share one form. Each node i has an activity xᵢ that evolves as

$$
\frac{dx_i}{dt} = M_0(x_i) + \sum_{j} A_{ij}\, M_1(x_i)\, M_2(x_j)
$$

M₀ is what a node does alone (decay, recovery), M₁ is how its own state gates incoming influence, and M₂ is what a neighbour sends. A is the adjacency matrix, and Sᵢ = Σⱼ Aᵢⱼ is the weighted degree.

I used three of the paper's models, one per propagation regime:

| Model | Equation | θ (theory) | Regime |
| --- | --- | --- | --- |
| ℝ₁, regulatory (a = 1, h = 1) | dx/dt = −x + Σ Aᵢⱼ xⱼ/(1 + xⱼ) | 0 | Distance-driven: arrival time counts hops |
| ℝ₂, regulatory (a = 0.4, h = 0.2) | dx/dt = −x⁰·⁴ + Σ Aᵢⱼ xⱼ⁰·²/(1 + xⱼ⁰·²) | 3/2 | Degree-driven: hubs are slow bottlenecks |
| 𝔼, SIS epidemic | dx/dt = −x + Σ Aᵢⱼ (1 − xᵢ) xⱼ | −1 | Composite: hubs are fast, low-degree nodes are slow |

The paper's central claim is that each node's response time scales as τᵢ \~ Sᵢ^θ, with θ fixed by the dynamics alone, and that arrival times across the network follow from summing these delays along paths.

## Theory: where θ comes from

θ follows from three steps: find each node's steady state, find how stiffly it is pulled back after a nudge, and read off how that stiffness scales with degree.

**Step 1, steady state.** Set dxᵢ/dt = 0 and use the mean-field approximation: in a random network every node's neighbours look statistically alike, so the neighbour sum is degree times a typical neighbour value m.

```latex
M_0(x_i) + M_1(x_i)\, S_i\, m = 0 \;\Longrightarrow\; x_i = x(S_i)
```

**Step 2, stiffness.** Hold the neighbours fixed after a small nudge and linearise. The node relaxes exponentially at a rate set by the diagonal Jacobian term, and the time to reach a fraction η of its final response is

```latex
|J_{ii}| = \left| M_0'(x_i) + M_1'(x_i)\, S_i\, m \right|, \qquad \tau_i = \frac{\ln\!\big(\tfrac{1}{1-\eta}\big)}{|J_{ii}|}
```

The size of the nudge cancels: it sets how far the node moves, not how fast.

**Step 3, scaling.** Substitute Step 1 into Step 2 and expand for large S. The leading power is θ.

- **SIS:** M₀′ = −B and M₁′ = −1 are constants, so |Jᵢᵢ| = B + Sᵢ m and τᵢ = 1.204/(B + Sᵢ m) at η = 0.7. At large S this is S⁻¹, so θ = −1. Hubs respond fast because every extra infection cuts the inflow from Sᵢ m infected neighbours.
- **ℝ₂:** Step 1 gives x = (S m)²·⁵, and the only restoring force is the slope of the decay term, 0.4 x⁻⁰·⁶ = 0.4 (S m)⁻¹·⁵. So τ ∝ S¹·⁵ and θ = 3/2. Hubs sit where the decay curve is nearly flat, so they barely get pulled back.
- **General regulatory model** with decay −xᵃ: θ = (1 − a)/a, which gives 0 for ℝ₁ and 3/2 for ℝ₂.

For SIS the power law is only asymptotic. The local log-log slope is

```latex
\frac{d\ln\tau}{d\ln S} = -\frac{S\,m}{B + S\,m}
```

which is about −0.41 at S = 1, −0.88 at S = 10 and −0.97 at S = 50. On a network whose degrees are mostly small, the measured slope sits well above −1.

**From one node to the network.** A signal travels as a relay, so delays add along the path. The paper's temporal distance (Eq. 5) sums Sₚ^θ over the nodes of a hop-shortest path, excluding the source, and minimises over such paths:

```latex
\mathcal{L}(j \to i) = \min_{\Pi(j \to i)} \sum_{p \in \Pi} S_p^{\theta}, \qquad T(j \to i) \propto \mathcal{L}(j \to i)
```

The claim is a proportionality: 𝓛 predicts relative arrival times, and one constant per system is fitted.

## Methods

Every experiment follows the paper's protocol: settle the system, clamp one source 10% higher, and time when each other node completes 70% of its response.

1. **Steady state.** Solve dxᵢ/dt = 0 on the actual network by fixed-point iteration (damped for ℝ₂). For SIS this is xᵢ = (Ax)ᵢ / (1 + (Ax)ᵢ).
2. **Perturbation.** Clamp the source at xⱼ = 1.1 xⱼ\* (α = 0.1) and integrate all other nodes from the old steady state.
3. **Final response.** Compute each node's new steady state directly, by fixed-point iteration with the source clamped, rather than trusting the last ODE sample. This rules out truncation errors.
4. **Arrival time.** With Δxᵢ(t) = xᵢ(t) − xᵢ\*, the normalised response is fᵢ(t) = Δxᵢ(t)/Δxᵢ(∞). T(j → i) is the first time fᵢ crosses η = 0.7, linearly interpolated between samples.
5. **Local response time.** τᵢ is the average of T(j → i) over the sampled sources j that are neighbours of i, as in the paper.

| Setting | Replication notebook | SIS and retrodiction notebook |
| --- | --- | --- |
| Networks | ER (giant component), BA scale-free with m = 2 | Same, plus SF1: the BA graph with weights U(0.1, 0.9) |
| Nodes N | 350 | 400 |
| Mean degree | 4 | 4 |
| Sources per network | 20 | 25 (all 400 for retrodiction) |
| Integrator | RK45 (ℝ₁, SIS); Radau to t = 6000 (ℝ₂) | RK45 to t = 15 |

The paper used N = 6,000. My networks are 15–17 times smaller, which matters most for SIS, where the power law only emerges at large degree.

## Replication: three regimes, three propagation patterns

The same networks give three completely different propagation patterns, and in each case the pattern matches the regime θ predicts.

**Figure 1** has one column per model (ℝ₁, ℝ₂, SIS) and one row per measurement.
<img width="1589" height="1289" alt="image" src="https://github.com/user-attachments/assets/2a854690-1182-4a92-a16f-cc69d9fc3080" />

*Figure 1. Columns: ℝ₁ (θ = 0), ℝ₂ (θ = 3/2), SIS (θ = −1). Rows: τᵢ vs Sᵢ; P(T) on the scale-free network; mean T vs hop distance. N = 350, 20 sources per model and network.*

**Row 1: response time τᵢ against degree Sᵢ (log-log).** Points are degree-binned averages; the black line is the theory slope θ.

- **ℝ₁:** flat at τ ≈ 1.2, which is ln(1/0.3) = 1.204 as derived. The apparent scatter is only about 2%; the y-axis spans 1.205 to 1.235.
- **ℝ₂:** rises with slope 3/2. A degree-1 node responds in about 1.7, a degree-50 hub in about 600. The mean-field formula τ = 3.01 (0.7 S)¹·⁵ gives 1.8 and 620, with no fitting.
- **SIS:** falls, but bends. The slope approaches −1 only at high degree, which is the crossover −S m/(B + S m) from the theory section.
- **ER and SF points fall on the same curve in every panel.** θ belongs to the dynamics, not the network. The network only decides which degrees exist.

**Row 2: distribution of arrival times P(T) on the scale-free network.**

- **ℝ₁:** a comb of sharp peaks. Each peak is one hop shell, spaced by about τ ≈ 1.2, with direct neighbours at the first spike. Arrival time literally counts hops.
- **ℝ₂:** times spread over three orders of magnitude and pile up just below 10³. That ceiling is set by the slowest hubs (τ ≈ 600). It is real, not a truncation: the largest T is about 700, the horizon is 6000, and no target is lost.
- **SIS:** a single smooth hump between 0 and 2.5. Fast hubs act as shortcuts, so hop shells blur together.

**Row 3: mean arrival time against hop distance dᵢⱼ.**

- **ℝ₁:** a straight line on both networks. Hop count explains everything.
- **ℝ₂:** ER rises gently. SF jumps at 2 hops and then falls with distance. More hops does not mean later arrival, because what matters is whether the path crosses a slow hub. Hop count is the wrong ruler.
- **SIS:** close to linear on both networks.

## Replication: does temporal distance predict arrival time?

T ∝ 𝓛 holds cleanly for both regulatory models on both networks, but not for SIS on the scale-free network.

**Figure 2** tests the paper's main claim directly. The x-axis is the temporal distance 𝓛(j → i); the y-axis is the measured arrival time divided by one fitted constant c per curve. Points are binned averages, and the dashed diagonal is a perfect prediction.
<img width="1590" height="449" alt="image" src="https://github.com/user-attachments/assets/913ec199-3fad-44ce-95ad-5b7e7af0ab94" />


*Figure 2. Collapse test, T/c against 𝓛(j → i) on log-log axes, binned. Dashed line: T = c𝓛.*

- **ℝ₁:** every point sits on the diagonal. With θ = 0, 𝓛 is just hop count, so this is a sanity check.
- **ℝ₂:** the paper's strongest case. Arrival times from 1 to about 700 land on the diagonal for both ER and SF. The non-monotonic curve in Figure 1, row 3 straightens into a line once distance is measured with 𝓛 instead of hops.
- **SIS:** ER mostly follows the line. On SF, targets with small 𝓛 sit far above the diagonal. Their paths run mostly through hubs, which S⁻¹ makes nearly free, yet the signal arrives much later than predicted.

Correlation of log T with each predictor, averaged over sources:

| Model | Network | Temporal distance 𝓛 | Hop count |
| --- | --- | --- | --- |
| ℝ₁ (θ = 0) | ER | 0.99 | 0.99 |
| ℝ₁ (θ = 0) | SF | 0.94 | 0.94 |
| ℝ₂ (θ = 3/2) | ER | 0.98 | 0.74 |
| ℝ₂ (θ = 3/2) | SF | **0.91** | **0.25** |
| SIS (θ = −1) | ER | 0.92 | 0.93 |
| SIS (θ = −1) | SF | **0.63** | **0.80** |

For ℝ₂ on SF, 𝓛 rescues a prediction that hop count gets almost entirely wrong. For SIS on SF the order reverses: the cruder metric wins.

## SIS in detail: where the metric breaks

The SIS failure on scale-free networks is not a finite-degree artefact: correcting the node weights to each node's exact stiffness closes only a small part of the gap.

**Figure 3** repeats the SIS measurements with a third network, SF1 (the same scale-free graph with edge weights drawn from U(0.1, 0.9)).
<img width="1590" height="470" alt="image" src="https://github.com/user-attachments/assets/6bca726c-25bc-457a-8356-f7d367d2fd6f" />

*Figure 3. SIS on ER, SF and SF1 (N = 400, 25 sources). Panel labels follow the paper's figure numbers; "Lᵢⱼ" in the middle and right panels denotes hop distance.*

- **Left, τᵢ against Sᵢ:** all three networks fall on one curve that bends away from the S⁻¹ line at low degree. SF1 extends below S = 1 because its weights are fractional, so it sits deepest in the crossover region.
- **Middle, P(T) split by hop shell:** the shells overlap heavily. A target 3 hops away can arrive before one 2 hops away, which is why the total distribution is a single hump rather than a comb.
- **Right, mean T against hop distance:** roughly linear at about 0.3 time units per hop on both networks. On average, hop count is a reasonable ruler for SIS.

**Figure 4** reproduces the paper's metric test and its wave snapshots.
<img width="1383" height="523" alt="image" src="https://github.com/user-attachments/assets/5dd89b13-a0c0-49a5-bb22-1e9b764434f1" />

*Figure 4. Left: 1/T vs 1/𝓛 for SIS, binned; solid line has slope 1. Middle and right: the scale-free network at t = 0.45 after perturbing a median-degree source, with radius = hop distance (middle) or 𝓛 (right); colour and size = fraction of final response.*

- **Left, 1/T against 1/𝓛:** the paper plots inverses, so short distances sit on the right. ER follows the diagonal. SF and SF1 fall below it at large 1/𝓛: targets the metric places close to the source actually respond much later than predicted. SF1 deviates most.
- **Middle and right, the network at t = 0.45:** each dot is a node, coloured and sized by how much of its final response it has completed. In the middle panel, radius is hop distance from the source; in the right panel, radius is 𝓛. Placing nodes by 𝓛 pulls the responded nodes toward the centre, so the picture looks more like a concentric wave. Visually this supports the paper, but the quantitative test on the left shows the fit is loose on SF.

**Exact-stiffness test.** For SIS the stiffness is known exactly with no mean-field step, |Jᵢᵢ| = 1 + (Ax\*)ᵢ. Replacing Sₚ⁻¹ with 1/|Jₚₚ| in 𝓛 tests whether the crossover explains the failure:

| Network | 𝓛 with S⁻¹ | 𝓛 with exact 1/J | Hop count | Steady-state x\* range |
| --- | --- | --- | --- | --- |
| ER | 0.92 | 0.95 | 0.93 | 0.34–0.89 |
| SF | 0.63 | 0.69 | 0.80 | 0.58–0.98 |

The exact weights help on both networks but leave 𝓛 well behind hop count on SF. The SF network is also much closer to saturation: its least-infected node already sits at x\* = 0.58.

## Retrodiction: finding the source from arrival times

For SIS on the scale-free network, arrival times identify the source almost perfectly, hop count recovers most of that, and the paper's temporal distance recovers essentially none of it.

**The problem.** One unknown node is perturbed. We see only the arrival times at a random set of observer nodes, possibly with noise, and we do not know exactly when the perturbation happened. Which node was the source?

**The estimator.** For each candidate source s, a predictor D(s, o) says how far in time each observer o is from s. If s is the true source, the observed times should line up as

```latex
t_o \approx t_0 + c\, D(s, o), \qquad c > 0
```

with unknown onset t₀ and unknown scale c. I fit both by least squares for every candidate, rank candidates by residual, and take the best as the estimate.

**Four predictors:**

- **hops:** D = hop distance, the classic source-localisation assumption.
- **𝓛 with S⁻¹:** the paper's temporal distance with θ = −1.
- **𝓛 with exact 1/J:** the same, with each node's exact SIS stiffness.
- **oracle:** D = the true simulated arrival times from s. This is the best any distance-based method could do, and it measures how much information the data contain.

**Protocol.** Every node of each 400-node network was simulated as a source. For each setting I ran 150 trials: random observers, a random true source that is never an observer, multiplicative timing noise of 0% or 20%, and a random onset shift of up to 5 time units.

**Figure 5** has four columns: ER and SF, each with 0% and 20% noise. The x-axis is the percentage of nodes observed (2, 5, 10, 20).

<img width="1790" height="790" alt="image" src="https://github.com/user-attachments/assets/05139d5e-f038-4ecf-8d98-f1ec0915764d" />


*Figure 5. SIS source retrodiction, 150 trials per point. Grey: hops. Green: 𝓛 with S⁻¹. Blue: 𝓛 with exact 1/J. Black dashed: oracle.*

- **Top row:** probability that the best guess is exactly the true source. Higher is better.
- **Bottom row:** average hop distance between the guess and the true source. Lower is better; random guessing gives about 3 hops.
- **ER:** every predictor improves with more observers, and they stay close together. Hops is slightly best, exact-stiffness 𝓛 next, S⁻¹ 𝓛 last.
- **SF:** the oracle reaches nearly 100% from 5–10% observed. Hops climbs steadily to 89%. Both 𝓛 estimators stay flat near 0%, missing by about 3 hops however many nodes are observed: no better than chance.

Top-1 accuracy with 20% of nodes observed:

| Network | Noise | Hops | 𝓛 with S⁻¹ | 𝓛 with exact 1/J | Oracle |
| --- | --- | --- | --- | --- | --- |
| ER | 0% | 0.93 | 0.83 | 0.92 | 1.00 |
| ER | 20% | 0.91 | 0.73 | 0.86 | 0.93 |
| SF | 0% | 0.89 | 0.01 | 0.03 | 1.00 |
| SF | 20% | 0.81 | 0.02 | 0.03 | 0.99 |

**Robustness check.** To rule out the linear fit as the cause, I refitted in log space with a known onset time, which matches the proportional form T ∝ 𝓛 exactly. On SF with 10% observed, 𝓛 rises only to 0.05 (S⁻¹) and 0.09 (exact 1/J), against 0.75 for hops and 0.99 for the oracle. The failure is in the metric, not the fitting.



3. Repeat at N = 6,000 and over several network realisations, and add SF1 and the paper's empirical networks.
4. Extend retrodiction to the degree-driven regime (ℝ₂), where 𝓛 should beat hop count, and map identifiability against θ.

*Code: two Jupyter notebooks, a replication notebook covering ℝ₁, ℝ₂ and SIS, and an SIS and retrodiction notebook.*
