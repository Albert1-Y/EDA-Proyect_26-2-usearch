#pragma once

#include <cstddef>
#include <cstdint>

#include "types.hpp"

namespace eda {

struct ResultadoSerializacion {
    char const* error = nullptr;

    explicit operator bool() const noexcept { return error == nullptr; }
};

struct ProgresoNulo {
    bool operator()(std::size_t, std::size_t) const noexcept { return true; }
};

struct CabeceraAlphaCG {
    char firma[8] = {'a', 'l', 'p', 'h', 'a', '-', 'c', 'g'};
    std::uint32_t version = 1;
    std::uint16_t bytes_slot = sizeof(slot_t); 
    std::uint16_t al_dia = 1;                  // 0: se guardo sin build(); se buscara por fuerza bruta
    std::uint64_t cantidad_nodos = 0;          // n
    std::uint64_t cantidad_aristas = 0;        // E
    std::uint64_t entrada = 0;
    std::uint64_t reservado = 0;
};
static_assert(sizeof(CabeceraAlphaCG) == 48, "la cabecera debe tener tamaño fijo y sin relleno");

} // namespace eda
