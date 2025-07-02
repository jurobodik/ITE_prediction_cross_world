"""
This script serves as a demo example to try out and possibly easy-to-modify for your dataset
First, we generate synthetic data using the `data_synthetic` function. It simulates realistic treatment-outcome dependencies using various copulas and marginal noise distributions.
Then, we apply the `CW_rho_intervals` function to estimate prediction intervals for Individual Treatment Effects (ITE).
The `CW_rho_intervals` function can be found in the `CW_rho_intervals_function.py` file

Function source: source(CW_rho_intervals_function.py)

This function is applied on the synthetic dataset and plotted the results below
"""

import numpy as np
import pandas as pd
from scipy.stats import norm, t, laplace, chi2, multivariate_normal
from copulas.multivariate import GaussianMultivariate
from copulas.bivariate import Clayton, Gumbel
from scipy.stats import multivariate_normal, t as t_dist



# Main synthetic data generation function
def data_synthetic(n=1000, d=2, rho=0, sigma_1=1, sigma_2=4, constant_propensity=False,
                   copula_type="gaussian", marginal="gaussian", random_function="perlin"):

    def random_function_1d(freq=0.05):
        if random_function == "perlin":
            from perlin_noise import PerlinNoise
            noise = PerlinNoise(octaves=2)
            s = np.linspace(-10, 10, 1001)
            f_noise = np.array([noise(x * freq) for x in s])
            amplitude = 30 / (np.max(f_noise) + 1)
            return f_noise * amplitude
        elif random_function == "sinusoidal":
            s = np.linspace(-10, 10, 1001)
            return np.sin(s * freq) * 10

    def evaluation_of_f_1d(f, X1):
        s_min, s_max = np.min(X1), np.max(X1)
        indices = ((X1 - s_min) * 1000 / (s_max - s_min)).astype(int)
        return f[indices]

    def random_function_2d(freq=0.001):
            x = np.linspace(-10, 10, 1001)
            y = np.linspace(-10, 10, 1001)
            X, Y = np.meshgrid(x, y)
            return np.sin(X * freq) * np.cos(Y * freq) * 10

    def evaluation_of_f_2d(f, X1, X2):
        x_min, x_max = np.min(X1), np.max(X1)
        y_min, y_max = np.min(X2), np.max(X2)
        x_idx = ((X1 - x_min) * 1000 / (x_max - x_min)).astype(int)
        y_idx = ((X2 - y_min) * 1000 / (y_max - y_min)).astype(int)
        return f[x_idx, y_idx]

    def generate_errors(n, rho, sigma_1, sigma_2, copula_type, marginal):
        marginal_dict = {
           "gaussian": norm.ppf,
           "t": lambda u: t.ppf(u, df=3),
            "laplace": laplace.ppf,
            "chisq": lambda u: chi2.ppf(u, df=3)
        }
        marginal_fn = marginal_dict[marginal]

        if copula_type == "gaussian":
            cov_matrix = np.array([[1, rho], [rho, 1]])
            u = norm.cdf(multivariate_normal.rvs(mean=[0, 0], cov=cov_matrix, size=n))
        elif copula_type == "t":
            cov_matrix = np.array([[1, rho], [rho, 1]])
            u = t_dist.cdf(t_dist.rvs(df=3, size=(n, 2)), df=3)
        elif copula_type in ["clayton", "gumbel"]:
            copula_dict = {"clayton": Clayton, "gumbel": Gumbel}
            copula = copula_dict[copula_type](theta=copula_dict[copula_type].tau_to_theta((2/np.pi)*np.arcsin(rho)))
            u = copula.sample(n)
        else:
            raise ValueError("Unsupported copula type.")

        eps1 = marginal_fn(u[:, 0]) * np.sqrt(sigma_1)
        eps2 = marginal_fn(u[:, 1]) * np.sqrt(sigma_2)

        return eps1, eps2

    if d == 1:
        X = np.random.uniform(-1, 1, n)
        CATE = 2+evaluation_of_f_1d(random_function_1d(), X)
        mu = 5 + 5 * X
    else:
        Sigma = np.full((d, d), 0.25)
        np.fill_diagonal(Sigma, 1)
        tilde_X = multivariate_normal.rvs(mean=np.zeros(d), cov=Sigma, size=n)
        X = norm.cdf(tilde_X)
        tau = random_function_2d()
        CATE = evaluation_of_f_2d(tau, X[:, 0], X[:, 1])
        beta = np.random.normal(size=d)
        mu = X @ beta

    epsilon1, epsilon2 = generate_errors(n, rho, sigma_1, sigma_2, copula_type, marginal)

    Y0 = mu + epsilon1
    Y1 = mu + CATE + epsilon2

    if not constant_propensity:
        propensity_score = (1 + np.abs(X if d == 1 else X[:, 0])) / 4
        treatment = np.random.binomial(1, 1 - propensity_score)
    else:
        treatment = np.random.binomial(1, 0.5, n)

    Y_obs = np.where(treatment == 1, Y1, Y0)

    data_dict = {'Y0': Y0, 'Y1': Y1, 'Y_obs': Y_obs, 'treatment': treatment}

    if d == 1:
        data_dict['X'] = X
    else:
        for i in range(d):
            data_dict[f'X{i+1}'] = X[:, i]

    return pd.DataFrame(data_dict)








import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split

rho = 0.1
#Generate synthetic data
df = data_synthetic(
    n=10000, d=1, rho=rho, sigma_1=1, sigma_2=4,
    constant_propensity=False, copula_type="gaussian",
    marginal="gaussian", random_function="perlin"
)

train_df, test_df = train_test_split(df, test_size=0.9) 

#Prepare inputs
X_train = pd.DataFrame(train_df.filter(like='X'))
Y_train = train_df['Y_obs'].values
treatment_train = train_df['treatment'].values

X_test =  pd.DataFrame(test_df.filter(like='X'))
Y0_test = test_df['Y0'].values
Y1_test = test_df['Y1'].values
true_ITE_test = Y1_test - Y0_test




#Run method on test points (but fit on train set!)
result = CW_rho_intervals(
    X=X_train,
    Y=Y_train,
    treatment=treatment_train,
    new_points=X_test.values,  # test X
    rho=rho,
    conformal='CQR',
    CQR_qr='qgam',     
    weighted_conformal=False,
    add_confidence_intervals=True
)

#Compute coverage: is true CATE in [lower, upper]?
lower = result['lower'].values
upper = result['upper'].values
pred_CATE = result['CATE'].values


coverage_rate = np.mean((true_ITE_test >= lower) & (true_ITE_test <= upper))
length = np.mean(upper - lower)
print(f"Coverage rate: {coverage_rate:.2f}")
print(f"Average length of intervals: {length:.2f}")



##Plot results if d==1
if X_train.shape[1] == 1:

  plt.figure(figsize=(12, 6))

  # Scatter observed Y_obs on test set (with treatment coloring)
  colors = test_df['treatment'].map({0: 'blue', 1: 'red'})
  plt.scatter(test_df['X'], test_df['Y_obs'], c=colors, alpha=0.5, s=20, label='Observed outcomes')

  # Sort for CATE line
  sorted_idx = X_test.values.flatten().argsort()
  X_sorted = X_test.values.flatten()[sorted_idx]
  cate_sorted = pred_CATE[sorted_idx]
  lower_sorted = lower[sorted_idx]
  upper_sorted = upper[sorted_idx]

  # Plot CATE and intervals
  plt.plot(X_sorted, cate_sorted, label='Estimated CATE', color='black', linewidth=2)
  plt.fill_between(X_sorted, lower_sorted, upper_sorted, color='green', alpha=0.3, label='90% Interval')

  # Plot true ITE points themselves
  plt.scatter(test_df['X'], true_ITE_test, color='green', alpha=0.1, s=10, label='True ITE points')

  # Add legend
  import matplotlib.patches as mpatches
  red_patch = mpatches.Patch(color='red', label='Treated (t=1)')
  blue_patch = mpatches.Patch(color='blue', label='Control (t=0)')

  plt.legend(handles=[
      red_patch, blue_patch, 
      plt.Line2D([0], [0], color='black', label='Estimated CATE'),
      mpatches.Patch(color='green', alpha=0.3, label='90% ITE Interval')
  ])

  # Labels and formatting
  plt.axhline(0, color='gray', linestyle='--')
  plt.xlabel("X")
  plt.ylabel("Y_obs / CATE")
  plt.title(f"ITE prediction intervals using ρ={rho}\nCoverage on test set: {coverage_rate:.3f}")
  plt.grid(True)
  plt.tight_layout()
  plt.show()

print(f"\nCoverage on test set: {coverage_rate:.3f}")




