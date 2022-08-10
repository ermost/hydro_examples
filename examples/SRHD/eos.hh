#pragma once

template <typename T = double> class SimpleGammaLaw {
  public:

  static T Gamma;

  template <typename Td>
  static inline T press_eps__rho_eps_th(T& eps_tot, Td const rho, Td const eps_th) {
    eps_tot = eps_th;
    return eps_tot * (Gamma - 1.);
  };

  template <typename Td>
  static inline T eps_th__eps_rho(Td const eps_tot, Td const rho) {
    return eps_tot;
  };

  template <typename Td> static inline T cs2__rho_eps_th(Td const rho, Td const eps_th) {
    auto const eps_tot = eps_th;
    auto press = eps_tot * (Gamma - 1.);

    auto rhoh = eps_tot*Gamma + rho;
    return Gamma * press / rhoh;
  };
};

template <typename T> T SimpleGammaLaw<T>::Gamma = 2.;
