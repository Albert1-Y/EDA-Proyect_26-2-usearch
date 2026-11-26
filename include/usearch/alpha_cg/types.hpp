#pragma once
// tipos basicos compartidos en los modulos

#include <cstddef>
#include <cstdint>

namespace eda {

using vector_key_t = std::uint64_t; // clave del usuario para un vector
using slot_t = std::uint32_t;       // indice interno del vector
using distance_t = float;

// vecino devuelto por el grafo
struct Neighbor {
    distance_t distancia;
    slot_t slot;
};

//primero distancia y si hay empate el slot
inline bool closer(Neighbor const& izquierdo, Neighbor const& derecho) noexcept {
    return izquierdo.distancia < derecho.distancia ||
           (izquierdo.distancia == derecho.distancia && izquierdo.slot < derecho.slot);
}

// contadores de una busqueda
struct ContadoresBusqueda {
    std::size_t distancias_calculadas = 0; 
    std::size_t nodos_expandidos = 0;
};

} // namespace eda
