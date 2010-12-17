from __future__ import print_function
from scipy.stats.distributions import norm as scipy_norm
from math import sqrt

norm = scipy_norm()

def pdf(x): return norm.pdf(x)
def cdf(x): return norm.cdf(x)
def icdf(x): return norm.ppf(x)    # inverse CDF

def Vwin(t, e):
  return pdf(t-e) / cdf(t-e)
def Wwin(t, e):
  return Vwin(t, e) * (Vwin(t, e) + t - e)

def Vdraw(t, e):
  return (pdf(-e-t) - pdf(e-t)) / (cdf(e-t) - cdf(-e-t))
def Wdraw(t, e):
  return Vdraw(t, e) ** 2 + ((e-t) * pdf(e-t) + (e+t) * pdf(e+t)) / (cdf(e-t) - cdf(-e-t))


class Gaussian(object):
  def __init__(self, mu=None, sigma=None, pi=None, tau=None):
    if pi is not None:
      self.pi = pi
      self.tau = tau
    elif mu is not None:
      self.pi = sigma ** -2
      self.tau = self.pi * mu
    else:
      self.pi = 0
      self.tau = 0

  def __repr__(self):
    return "N(pi={0.pi},tau={0.tau})".format(self)

  def __str__(self):
    if self.pi == 0.0:
      return "N(mu=0,sigma=inf)"
    else:
      sigma = sqrt(1/self.pi)
      mu = self.tau / self.pi
      return "N(mu={0:.3f},sigma={1:.3f})".format(mu, sigma)

  def MuSigma(self):
    if self.pi == 0.0:
      return 0, float("inf")
    else:
      return self.tau / self.pi, sqrt(1/self.pi)


  def __mul__(self, other):
    return Gaussian(pi=self.pi+other.pi, tau=self.tau+other.tau)

  def __div__(self, other):
    return Gaussian(pi=self.pi-other.pi, tau=self.tau-other.tau)


class Variable(object):
  def __init__(self):
    self.value = Gaussian()
    self.factors = {}

  def AttachFactor(self, factor):
    self.factors[factor] = Gaussian()

  def UpdateMessage(self, factor, message):
    old_message = self.factors[factor]
    self.value = self.value / old_message * message
    self.factors[factor] = message

  def UpdateValue(self, factor, value):
    old_message = self.factors[factor]
    self.factors[factor] = value * old_message / self.value
    self.value = value

  def GetMessage(self, factor):
    return self.factors[factor]


class Factor(object):
  def __init__(self, variables):
    self.variables = variables
    self.value = Gaussian()
    for v in variables:
      v.AttachFactor(self)


class PriorFactor(Factor):
  def __init__(self, variable, param):
    super(PriorFactor, self).__init__([variable])
    self.param = param

  def Start(self):
    self.variables[0].UpdateValue(self, self.param)


class LikelihoodFactor(Factor):
  def __init__(self, mean_variable, value_variable, variance):
    super(LikelihoodFactor, self).__init__([mean_variable, value_variable])
    self.mean = mean_variable
    self.value = value_variable
    self.variance = variance

  def UpdateValue(self):
    y = self.mean.value
    f_y = self.mean.GetMessage(self)
    a = 1.0 / (1.0 + self.variance * (y.pi - f_y.pi))
    self.value.UpdateMessage(self, Gaussian(pi=a*(y.pi - f_y.pi),
                                            tau=a*(y.tau - f_y.tau)))

  def UpdateMean(self):
    y = self.value.value
    f_y = self.value.GetMessage(self)
    a = 1.0 / (1.0 + self.variance * (y.pi - f_y.pi))
    self.mean.UpdateMessage(self, Gaussian(pi=a*(y.pi - f_y.pi),
                                           tau=a*(y.tau - f_y.tau)))


class SumFactor(Factor):
  def __init__(self, sum_variable, terms_variables, coeffs):
    assert len(terms_variables) == len(coeffs)
    self.sum = sum_variable
    self.terms = terms_variables
    self.coeffs = coeffs
    super(SumFactor, self).__init__([sum_variable] + terms_variables)

  def _InternalUpdate(self, var, y, fy, a):
    new_pi = 1.0 / (sum(a[j]**2 / (y[j].pi - fy[j].pi) for j in range(len(a))))
    new_tau = new_pi * sum(a[j] *
                           (y[j].tau - fy[j].tau) / (y[j].pi - fy[j].pi)
                           for j in range(len(a)))
    var.UpdateMessage(self, Gaussian(pi=new_pi, tau=new_tau))

  def UpdateSum(self):
    y = [t.value for t in self.terms]
    fy = [t.GetMessage(self) for t in self.terms]
    a = self.coeffs
    self._InternalUpdate(self.sum, y, fy, a)

  def UpdateTerm(self, index):
    b = self.coeffs
    a = [-b[i] / b[index] for i in range(len(b)) if i != index]
    a.insert(index, 1.0 / b[index])

    v = self.terms[:]
    v[index] = self.sum
    y = [i.value for i in v]
    fy = [i.GetMessage(self) for i in v]
    self._InternalUpdate(self.terms[index], y, fy, a)


class TruncateFactor(Factor):
  def __init__(self, variable, V, W, epsilon):
    super(TruncateFactor, self).__init__([variable])
    self.var = variable
    self.V = V
    self.W = W
    self.epsilon = epsilon

  def Update(self):
    x = self.var.value
    fx = self.var.GetMessage(self)

    c = x.pi - fx.pi
    d = x.tau - fx.tau
    sqrt_c = sqrt(c)
    args = (d / sqrt_c, self.epsilon * sqrt_c)
    V = self.V(*args)
    W = self.W(*args)
    new_val = Gaussian(pi=c / (1.0 - W), tau=(d + sqrt_c * V) / (1.0 - W))
    self.var.UpdateValue(self, new_val)


def DrawProbability(epsilon, beta, total_players=2):
  return 2 * cdf(epsilon / (sqrt(total_players) * beta)) - 1

def DrawMargin(p, beta, total_players=2):
  return icdf((p+1.0)/2) * sqrt(total_players) * beta


INITIAL_MU = 25.0
INITIAL_SIGMA = INITIAL_MU / 3.0
BETA = INITIAL_SIGMA / 2.0
EPSILON = DrawMargin(0.05, BETA)   # 5% draw probability


def AdjustPlayers(players):
  players = players[:]
  players.sort(key=lambda p: p.rank)

  ss = [Variable() for p in players]
  ps = [Variable() for p in players]
  ts = [Variable() for p in players]
  ds = [Variable() for p in players[:-1]]

  skill = [PriorFactor(s, Gaussian(mu=pl.skill[0], sigma=pl.skill[1]))
           for (s, pl) in zip(ss, players)]
  skill_to_perf = [LikelihoodFactor(s, p, BETA**2)
                   for (s, p) in zip(ss, ps)]
  perf_to_team = [SumFactor(t, [p], [1])
                  for (p, t) in zip(ps, ts)]
  team_diff = [SumFactor(d, [t1, t2], [+1, -1])
               for (d, t1, t2) in zip(ds, ts[:-1], ts[1:])]
  trunc = [TruncateFactor(d, Vdraw if pl1.rank == pl2.rank else Vwin,
                          Wdraw if pl1.rank == pl2.rank else Wwin, EPSILON)
           for (d, pl1, pl2) in zip(ds, players[:-1], players[1:])]

  for f in skill:
    f.Start()
  for f in skill_to_perf:
    f.UpdateValue()
  for f in perf_to_team:
    f.UpdateSum()

  for i in range(5):
    for f in team_diff:
      f.UpdateSum()
    for f in trunc:
      f.Update()
    for f in team_diff:
      f.UpdateTerm(0)
      f.UpdateTerm(1)

  for f in perf_to_team:
    f.UpdateTerm(0)
  for f in skill_to_perf:
    f.UpdateMean()

  for s, pl in zip(ss, players):
    pl.skill = s.value.MuSigma()


__all__ = ["AdjustPlayers", "INITIAL_MU", "INITIAL_SIGMA"]
