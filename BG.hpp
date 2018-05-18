#ifndef _NUMCHY_BG_HPP
#define _NUMCHY_BG_HPP 11

#include <array>
#include <vector>
#include "lorentz.hpp"
#include "arra.hpp"

using athl::L4v;

//color ordered amplitude using BG recurrance relations
std::complex<double> amp_CO_BG(const std::vector<L4v<double>>& p, std::vector<int> h);

using athl::arra;
using athl::acca;

acca amp_CO_BG(const std::vector<L4v<arra>>& p, std::vector<int> h);

#endif
