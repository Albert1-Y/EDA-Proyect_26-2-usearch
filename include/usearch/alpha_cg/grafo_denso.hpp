#pragma once

// Fase 1 del grafo alpha-CG
//  grafo_denso_exacto      n^2 distancias.
//  grafo_denso_nn_descent  n^1.14

#include <algorithm>
#include <cstddef>
#include <random>
#include <vector>

#include "types.hpp"

namespace eda {

struct ParametrosNNDescent {
    float muestreo = 0.5f;                 
    float umbral_convergencia = 0.001f;
    std::size_t iteraciones_maximas = 12;
    unsigned semilla = 42;  
    std::size_t usar_fuerza_bruta_hasta = 4;
};

struct EstadisticasGrafoDenso {
    bool con_nn_descent = false;
    std::size_t iteraciones = 0;
    std::size_t distancias_calculadas = 0;
    std::vector<std::size_t> cambios_por_iteracion; 

// fuerza bruta
template <class Distancias>
std::vector<std::vector<slot_t>> grafo_denso_exacto(std::size_t cantidad_nodos, std::size_t k,
                                                    Distancias const& distancias,
                                                    EstadisticasGrafoDenso* estadisticas = nullptr)
{
    std::vector<std::vector<slot_t>> grafo(cantidad_nodos);
    std::vector<Neighbor> todos;
    for (slot_t p = 0; p != cantidad_nodos; ++p) {
        todos.clear();
        for (slot_t otro = 0; otro != cantidad_nodos; ++otro)
            if (otro != p)
                todos.push_back( {distancias.entre(p, otro), otro} );
        std::size_t const cuantos = (std::min)(k, todos.size());
        std::partial_sort(todos.begin(), todos.begin() + cuantos, todos.end(), closer);
        for (std::size_t i = 0; i != cuantos; ++i)
            grafo[p].push_back(todos[i].slot);
    }
    if (estadisticas) {
        *estadisticas = {};
        estadisticas->distancias_calculadas = cantidad_nodos * (cantidad_nodos ? cantidad_nodos - 1 : 0);
    }
    return grafo;
}


// NN-Descent el vecino de mi vecino probablemente sea mi vecino
namespace detalle_nn_descent {

struct VecinoMarcado {
    distance_t distancia;
    slot_t slot;
    bool nuevo; // local join
};

using ListaKnn = std::vector<VecinoMarcado>;

// intenta insertar un vecino
inline bool intentar_insertar(ListaKnn& lista, std::size_t k, distance_t distancia, slot_t slot)
{
    Neighbor const candidato{distancia, slot};
    bool const lista_llena = lista.size() >= k;

    if ( lista_llena && !closer( candidato, Neighbor{lista.back().distancia, lista.back().slot} ) )
        return false; // no mejora a nadie
    
    for (VecinoMarcado const& existente : lista)
        if (existente.slot == slot)
            return false;
    
    // insertar en orden
    auto posicion = lista.begin();
    while (posicion != lista.end() &&
        posicion->distancia < candidato.distancia)
    {
        ++posicion;
    }
    lista.insert(posicion, {distancia, slot, true});

    if (lista.size() > k)
        lista.pop_back(); // sale el mas lejano
    return true;
}

inline void muestrear(std::vector<slot_t>& lista, std::size_t maximo, std::mt19937& azar) {
    if (lista.size() <= maximo)
        return;
    std::shuffle(lista.begin(), lista.end(), azar);
    lista.resize(maximo);
}

inline void quitar_repetidos(std::vector<slot_t>& lista) {
    std::sort(lista.begin(), lista.end());
    lista.erase( std::unique( lista.begin(), lista.end() ), lista.end() );
}

} // namespace detalle_nn_descent

template <class Distancias>
std::vector<std::vector<slot_t>> grafo_denso_nn_descent(std::size_t cantidad_nodos, std::size_t k,
                                                        Distancias const& distancias,
                                                        ParametrosNNDescent const& parametros = {},
                                                        EstadisticasGrafoDenso* estadisticas = nullptr)
{

    if (cantidad_nodos <= parametros.usar_fuerza_bruta_hasta * k + 1)
        return grafo_denso_exacto(cantidad_nodos, k, distancias, estadisticas);

    EstadisticasGrafoDenso registro;
    registro.con_nn_descent = true;

    std::mt19937 azar(parametros.semilla);
    std::vector<detalle_nn_descent::ListaKnn> listas(cantidad_nodos);

    auto distancia = [&](slot_t a, slot_t b) {
        ++registro.distancias_calculadas;
        return distancias.entre(a, b);
    };

    // paso 1: K vecinos al azar para cada nodo
    for (slot_t v = 0; v != cantidad_nodos; ++v)
    {   
        while (listas[v].size() < k)
        {
            slot_t const otro = (slot_t)(azar() % cantidad_nodos);
            if (otro != v)
                detalle_nn_descent::intentar_insertar(listas[v], k, distancia(v, otro), otro);
        }
    }

    std::size_t const muestra = (std::max)( std::size_t{1}, (std::size_t)( parametros.muestreo * (float)(k) ) );
    
    double const minimo_de_cambios = parametros.umbral_convergencia * (double)(cantidad_nodos) * (double)(k);

    // paso 2 local join
    for (std::size_t ronda = 0; ronda != parametros.iteraciones_maximas; ++ronda) {
        
        // separar nuevos y viejos
        std::vector< std::vector<slot_t> > nuevos(cantidad_nodos), viejos(cantidad_nodos);

        for (slot_t v = 0; v != cantidad_nodos; ++v)
        {            
            std::vector<std::size_t> posiciones_nuevas;
            for (std::size_t i = 0; i != listas[v].size(); ++i)
            {
                if (listas[v][i].nuevo)
                    posiciones_nuevas.push_back(i);
                else
                    viejos[v].push_back(listas[v][i].slot);
            }

            std::shuffle(posiciones_nuevas.begin(), posiciones_nuevas.end(), azar);
            if (posiciones_nuevas.size() > muestra)
                posiciones_nuevas.resize(muestra);
            
            for (std::size_t i : posiciones_nuevas) {
                nuevos[v].push_back(listas[v][i].slot);
                listas[v][i].nuevo = false;
            }
        }

        // vecinos inversos
        std::vector<std::vector<slot_t>> nuevos_inversos(cantidad_nodos), viejos_inversos(cantidad_nodos);

        for (slot_t v = 0; v != cantidad_nodos; ++v)
        {
            for (slot_t u : nuevos[v])
                nuevos_inversos[u].push_back(v);
            for (slot_t u : viejos[v])
                viejos_inversos[u].push_back(v);
        }

        for (slot_t v = 0; v != cantidad_nodos; ++v)
        {
            detalle_nn_descent::muestrear(nuevos_inversos[v], muestra, azar);
            detalle_nn_descent::muestrear(viejos_inversos[v], muestra, azar);
            
            nuevos[v].insert(nuevos[v].end(), nuevos_inversos[v].begin(), nuevos_inversos[v].end());
            viejos[v].insert(viejos[v].end(), viejos_inversos[v].begin(), viejos_inversos[v].end());
            
            detalle_nn_descent::quitar_repetidos(nuevos[v]);
            detalle_nn_descent::quitar_repetidos(viejos[v]);
        }

        // local join
        std::size_t cambios = 0;

        auto presentar = [&](slot_t a, slot_t b) {
            if (a == b)
                return;
            distance_t const d = distancia(a, b);
            if (detalle_nn_descent::intentar_insertar(listas[a], k, d, b))
                ++cambios;
            if (detalle_nn_descent::intentar_insertar(listas[b], k, d, a))
                ++cambios;
        };

        for (slot_t v = 0; v != cantidad_nodos; ++v) {
            auto const& n = nuevos[v];
            for (std::size_t i = 0; i != n.size(); ++i)
                for (std::size_t j = i + 1; j != n.size(); ++j)
                    presentar(n[i], n[j]);       // nuevo - nuevo
            for (slot_t a : n)
                for (slot_t b : viejos[v])
                    presentar(a, b);             // nuevo - viejo
            // viejo - viejo ya se presentaron siendo nuevos
        }

        registro.iteraciones = ronda + 1;
        registro.cambios_por_iteracion.push_back(cambios);
        if ( (double)(cambios) < minimo_de_cambios )
            break;
    }

    // paso 3 devolver solo los slots
    std::vector<std::vector<slot_t>> grafo(cantidad_nodos);
    
    for (slot_t v = 0; v != cantidad_nodos; ++v)
    {
        for (detalle_nn_descent::VecinoMarcado const& vecino : listas[v])
            grafo[v].push_back(vecino.slot);
    }

    if (estadisticas)
        *estadisticas = std::move(registro);
    return grafo;
}

} // namespace eda
