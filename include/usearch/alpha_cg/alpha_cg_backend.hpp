#pragma once

// fase 2 del grafo alpha-CG

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <limits>
#include <vector>

#include "grafo_denso.hpp"
#include "serializacion.hpp"
#include "types.hpp"

namespace eda {

struct ParametrosAlphaCG {
    std::size_t vecinos_densos = 64;          // Paper: 200
    bool usar_nn_descent = true;
    ParametrosNNDescent nn_descent;
    std::size_t max_vecinos = 32;             // M
    std::size_t ancho_beam_construccion = 64; // L
    std::size_t max_candidatos = 500;         // C candidatos a podar
    float alpha_inicial = 0.9f;               // poda adaptativa
    float alpha_maximo = 1.6f;
    float paso_alpha = 0.05f;
    float tau = 0.0f;
    std::size_t ancho_beam_busqueda = 64;
    std::size_t muestra_para_entrada = 256;   // nodos que se comparan para elegir la entrada (fase 2)
};

// poda alpha-CG
template <class DistanciaEntre>
std::vector<Neighbor> podar(std::vector<Neighbor> const& candidatos, float alpha, float tau, DistanciaEntre& distancia_entre)
{
    std::vector<Neighbor> elegidos;
    for (Neighbor const& u : candidatos) {
        bool cubierto = false;
        for (Neighbor const& v : elegidos) {
            float const distancia_p_u = u.distancia;
            float const distancia_u_v = distancia_entre(u.slot, v.slot);
            if (distancia_p_u > alpha * distancia_u_v + (alpha + 1.0f) * tau) {
                cubierto = true; // se llega a u pasando por v
                break;
            }
        }
        if (!cubierto)
            elegidos.push_back(u);
    }
    return elegidos;
}

// oda adaptativa 
template <class DistanciaEntre>
std::vector<Neighbor> poda_adaptativa(std::vector<Neighbor> const& candidatos, ParametrosAlphaCG const& parametros,
                                      DistanciaEntre& distancia_entre)
{
    std::size_t const minimo_deseado = parametros.max_vecinos / 2;
    std::vector<Neighbor> elegidos;
    for (float alpha = parametros.alpha_inicial; alpha <= parametros.alpha_maximo + 1e-6f; alpha += parametros.paso_alpha)
    {
        elegidos = podar(candidatos, alpha, parametros.tau, distancia_entre);
        if (elegidos.size() >= minimo_deseado)
            break;
    }
    if (elegidos.size() > parametros.max_vecinos)
        elegidos.resize(parametros.max_vecinos); // quedan los mas cercanos
    return elegidos;
}


// busqueda beam
template <class DistanciaA>
std::vector<Neighbor> recorrer_beam(slot_t entrada, std::size_t ancho,
                                    std::vector<std::vector<slot_t>> const& grafo, DistanciaA&& distancia_a,
                                    ContadoresBusqueda* contadores = nullptr)
{
    struct EnBeam {
        Neighbor nodo;
        bool expandido; // true si ya recorrimos sus vecinos
    };

    std::vector<EnBeam> beam;
    std::vector<char> ya_visto( grafo.size(), 0 );
    std::vector<Neighbor> visitados;

    auto visitar = [&](slot_t slot)
    {
        ya_visto[slot] = 1;

        float distancia = distancia_a(slot);
        if (contadores)
            ++contadores->distancias_calculadas;
        
        Neighbor nodo{distancia, slot};
        visitados.push_back(nodo);

        bool beam_lleno = beam.size() >= ancho;
        if (beam_lleno)
        {
            Neighbor nodo_mas_lejano = beam.back().nodo;
            bool es_mejor = closer(nodo, nodo_mas_lejano);
            if (!es_mejor) return; // no entra al beam
        }

        auto posicion = beam.begin();
        while (posicion != beam.end())
        {
            bool nodo_actual_es_mas_lejano = closer(nodo, posicion->nodo);
            if (nodo_actual_es_mas_lejano) break;
            ++posicion;
        }

        beam.insert(posicion, EnBeam{nodo, false});
        if (beam.size() > ancho)
            beam.pop_back();
    };

    visitar(entrada);
    while (true)
    {
        auto siguiente = beam.begin();
        while (siguiente != beam.end() && siguiente->expandido) ++siguiente;

        if (siguiente == beam.end()) break;

        siguiente->expandido = true;
        if (contadores)
            ++contadores->nodos_expandidos;

        slot_t actual = siguiente->nodo.slot;

        for (slot_t vecino : grafo[actual])
        {
            if (ya_visto[vecino] == 0)
                visitar(vecino);
        }
    }
    return visitados;
}


class AlphaCGBackend {
  public:
    struct Estadisticas {
        EstadisticasGrafoDenso grafo_denso;
        double grado_promedio_tras_poda = 0;
        double grado_promedio_final = 0;
        std::size_t grado_maximo_final = 0;
        std::size_t aristas_inversas_agregadas = 0;
        std::size_t nodos_repodados = 0;
        std::size_t nodos_reconectados = 0;
    };

    AlphaCGBackend() = default;
    explicit AlphaCGBackend(ParametrosAlphaCG parametros) : parametros_(parametros) {}

    void reserve(std::size_t cantidad_slots) { vecinos_.reserve(cantidad_slots); }

    void add(slot_t slot) {
        assert(slot == cantidad_slots_ && "los slots nuevos deben llegar en orden");
        (void)slot;
        ++cantidad_slots_;
        grafo_al_dia_ = false;
    }

    void update(slot_t slot) {
        assert(slot < cantidad_slots_);
        (void)slot;
        grafo_al_dia_ = false; // el grafo ya no es valido
    }

    // construye el grafo completo
    template <class Distancias>
    void build(Distancias const& distancias)
    {
        std::size_t const cantidad_nodos = cantidad_slots_;
        vecinos_.assign(cantidad_nodos, {});
        estadisticas_ = {};
        entrada_ = 0;

        if (cantidad_nodos == 0) {
            grafo_al_dia_ = true;
            return;
        }
        
        // L2Sq -> L2, coseno -> cuerda
        auto distancia_entre = [&](slot_t a, slot_t b) {
            return distancias.real(distancias.entre(a, b));
        };

        std::vector<std::vector<slot_t>> const grafo_denso = fase1_grafo_denso(distancias);

        entrada_ = fase2_elegir_entrada(distancias);

        fase3_podar_cada_nodo(grafo_denso, distancias, distancia_entre);
        estadisticas_.grado_promedio_tras_poda = grado_promedio();
        fase4_agregar_inversas(distancia_entre);
        fase5_asegurar_conectividad(distancias);

        estadisticas_.grado_promedio_final = grado_promedio();
        for (auto const& lista : vecinos_)
            estadisticas_.grado_maximo_final = (std::max)(estadisticas_.grado_maximo_final, lista.size());
        grafo_al_dia_ = true;
    }

    template <class Distancias, class Allow>
    std::vector<Neighbor> search(typename Distancias::consulta_t consulta, std::size_t k,
                                 std::size_t ancho_beam, Distancias const& distancias,
                                 Allow& permitir, ContadoresBusqueda* contadores = nullptr) const
    {
        if (contadores)
            *contadores = {}; // cuenta solo esta busqueda
        std::vector<Neighbor> encontrados;
        if (!k || cantidad_slots_ == 0)
            return encontrados;

        if (!grafo_al_dia_) {
            // si el grafo no esta al dia -> fuerza bruta
            for (slot_t slot = 0; slot != cantidad_slots_; ++slot)
                if (permitir(slot))
                    encontrados.push_back( { distancias.desde(consulta, slot), slot } );
            if (contadores)
                contadores->distancias_calculadas = encontrados.size();
        }
        else
        {
            std::size_t const ancho = (std::max)(k, ancho_beam ? ancho_beam : parametros_.ancho_beam_busqueda);
            auto const distancia_a_consulta = [&](slot_t slot) { return distancias.desde(consulta, slot); };

            for ( Neighbor const& nodo : recorrer_beam(entrada_, ancho, vecinos_, distancia_a_consulta, contadores) )
                if (permitir(nodo.slot))
                    encontrados.push_back(nodo);
        }
        std::size_t const cuantos = (std::min)(k, encontrados.size());
        std::partial_sort(encontrados.begin(), encontrados.begin() + cuantos, encontrados.end(), closer);
        encontrados.resize(cuantos);
        return encontrados;
    }

    std::size_t size() const noexcept { return cantidad_slots_; }

    void clear() noexcept {
        cantidad_slots_ = 0;
        vecinos_.clear();
        entrada_ = 0;
        grafo_al_dia_ = true;
        estadisticas_ = {};
    }

    // serializacion

    std::size_t serialized_length() const noexcept {
        return sizeof(CabeceraAlphaCG) + (cantidad_slots_ + std::size_t{1}) * sizeof(std::uint64_t) +
               cantidad_aristas() * sizeof(slot_t);
    }

    template <class Escribir, class Progreso = ProgresoNulo>
    ResultadoSerializacion save_to_stream(Escribir&& escribir, Progreso&& progreso = {}) const {
        std::size_t const n = cantidad_slots_;

        CabeceraAlphaCG cabecera;
        cabecera.al_dia = grafo_al_dia_ ? 1 : 0;
        cabecera.cantidad_nodos = n;
        cabecera.cantidad_aristas = cantidad_aristas();
        cabecera.entrada = entrada_;
        if (!escribir(&cabecera, sizeof(cabecera)))
            return {"no se pudo escribir la cabecera"};

        std::uint64_t inicio = 0;
        if (!escribir(&inicio, sizeof(inicio)))
            return {"no se pudo escribir inicio[]"};
        for (slot_t v = 0; v != n; ++v) {
            inicio += lista_de(v).size();
            if (!escribir(&inicio, sizeof(inicio)))
                return {"no se pudo escribir inicio[]"};
        }

        for (slot_t v = 0; v != n; ++v) {
            std::vector<slot_t> const& lista = lista_de(v);
            if (!lista.empty() && !escribir(lista.data(), lista.size() * sizeof(slot_t)))
                return {"no se pudo escribir vecinos[]"};
            if (!progreso(std::size_t{v} + 1, n))
                return {"cancelado por el usuario"};
        }
        return {};
    }

    template <class Leer, class Progreso = ProgresoNulo>
    ResultadoSerializacion load_from_stream(Leer&& leer, Progreso&& progreso = {}) {
        clear();
        auto fallo = [&](char const* mensaje) {
            clear();
            return ResultadoSerializacion{mensaje};
        };

        CabeceraAlphaCG cabecera;
        if (!leer(&cabecera, sizeof(cabecera)))
            return fallo("no se pudo leer la cabecera");
        if (std::memcmp(cabecera.firma, CabeceraAlphaCG{}.firma, sizeof(cabecera.firma)) != 0)
            return fallo("no es un grafo alpha-CG");
        if (cabecera.version != CabeceraAlphaCG{}.version)
            return fallo("version de archivo no soportada");
        if (cabecera.bytes_slot != sizeof(slot_t))
            return fallo("el archivo usa otro tamaño de slot");

        std::uint64_t const n = cabecera.cantidad_nodos;
        std::uint64_t const aristas = cabecera.cantidad_aristas;
        if (n > (std::numeric_limits<slot_t>::max)())
            return fallo("demasiados nodos para slot_t");
        if (n == 0 ? aristas != 0 : cabecera.entrada >= n)
            return fallo("archivo corrupto: entrada o aristas fuera de rango");

        // inicio[]: debe empezar en 0, no decrecer y terminar en E
        std::vector<std::uint64_t> inicio(static_cast<std::size_t>(n) + 1);
        if (!leer(inicio.data(), inicio.size() * sizeof(std::uint64_t)))
            return fallo("no se pudo leer inicio[]");
        if (inicio.front() != 0 || inicio.back() != aristas)
            return fallo("archivo corrupto: inicio[] no cuadra con E");
        for (std::size_t v = 0; v != n; ++v)
            if (inicio[v + 1] < inicio[v])
                return fallo("archivo corrupto: inicio[] decrece");

        std::vector<std::vector<slot_t>> vecinos(static_cast<std::size_t>(n));
        for (std::size_t v = 0; v != n; ++v) {
            std::vector<slot_t>& lista = vecinos[v];
            lista.resize(static_cast<std::size_t>(inicio[v + 1] - inicio[v]));
            if (!lista.empty() && !leer(lista.data(), lista.size() * sizeof(slot_t)))
                return fallo("no se pudo leer vecinos[]");
            for (slot_t vecino : lista)
                if (vecino >= n)
                    return fallo("archivo corrupto: vecino fuera de rango");
            if (!progreso(v + 1, static_cast<std::size_t>(n)))
                return fallo("cancelado por el usuario");
        }

        cantidad_slots_ = static_cast<slot_t>(n);
        vecinos_ = std::move(vecinos);
        entrada_ = static_cast<slot_t>(cabecera.entrada);
        grafo_al_dia_ = cabecera.al_dia != 0 || n == 0;
        return {};
    }

    ResultadoSerializacion view(void const*, std::size_t /*longitud*/, std::size_t /*offset*/ = 0) {
        return {"view todavia no esta soportado: usa load_from_stream"};
    }

    // para estudiar el grafo

    bool grafo_al_dia() const noexcept { return grafo_al_dia_; }
    slot_t entrada() const noexcept { return entrada_; }
    std::vector<slot_t> const& vecinos_de(slot_t slot) const { return vecinos_[slot]; }
    ParametrosAlphaCG const& parametros() const noexcept { return parametros_; }
    Estadisticas const& estadisticas() const noexcept { return estadisticas_; }

  private:
    
    template <class Distancias>
    std::vector<std::vector<slot_t>> fase1_grafo_denso(Distancias const& distancias)
    {
        if (parametros_.usar_nn_descent)
            return grafo_denso_nn_descent(cantidad_slots_, parametros_.vecinos_densos, distancias,
                                          parametros_.nn_descent, &estadisticas_.grafo_denso);
        return grafo_denso_exacto(cantidad_slots_, parametros_.vecinos_densos, distancias,
                                  &estadisticas_.grafo_denso);
    }

    // fase 2 -> entrada
    template <class Distancias>
    slot_t fase2_elegir_entrada(Distancias const& distancias) const
    {
        std::size_t const tamano_muestra = (std::min)((std::size_t)cantidad_slots_, parametros_.muestra_para_entrada);
        std::vector<slot_t> muestra;
        for (std::size_t i = 0; i != tamano_muestra; ++i)
            muestra.push_back( (slot_t)(i * cantidad_slots_ / tamano_muestra) ); // repartidos parejo

        Neighbor mejor{0, muestra[0]};
        for (std::size_t i = 0; i != muestra.size(); ++i)
        {
            Neighbor candidato{0, muestra[i]};
            for (slot_t otro : muestra)
                candidato.distancia += distancias.entre(muestra[i], otro);
            if (i == 0 || closer(candidato, mejor))
                mejor = candidato;
        }
        return mejor.slot;
    }

    // fase 3 -> podar cada nodo
    template <class Distancias, class DistanciaEntre>
    void fase3_podar_cada_nodo(std::vector<std::vector<slot_t>> const& grafo_denso,
                               Distancias const& distancias, DistanciaEntre& distancia_entre)
    {
        for (slot_t p = 0; p != cantidad_slots_; ++p)
        {
            auto const distancia_a_p = [&](slot_t slot) { return distancias.entre(p, slot); };
            std::vector<Neighbor> candidatos = recorrer_beam(entrada_, parametros_.ancho_beam_construccion,
                                                             grafo_denso, distancia_a_p);

            for (slot_t vecino : grafo_denso[p])
                candidatos.push_back( {distancias.entre(p, vecino), vecino} );

            // ordenamos por slot para eliminar duplicados
            std::sort(candidatos.begin(), candidatos.end(), [](Neighbor const& a, Neighbor const& b) { return a.slot < b.slot; });
            // eliminamos repetidos
            candidatos.erase(std::unique(candidatos.begin(), candidatos.end(),
                                        [](Neighbor const& a, Neighbor const& b) { return a.slot == b.slot; }),
                             candidatos.end());
            // eliminamos p
            candidatos.erase(std::remove_if(candidatos.begin(), candidatos.end(),
                                            [p](Neighbor const& c) { return c.slot == p; }),
                             candidatos.end());
            
            for (Neighbor& c : candidatos)
                c.distancia = distancias.real(c.distancia); // a distancia metrica (L2Sq -> L2)
            
            std::sort(candidatos.begin(), candidatos.end(), closer);

            if (candidatos.size() > parametros_.max_candidatos)
                candidatos.resize(parametros_.max_candidatos);

            // podar los candidatos elegidos
            for (Neighbor const& elegido : poda_adaptativa(candidatos, parametros_, distancia_entre))
                vecinos_[p].push_back(elegido.slot);
        }
    }

    // fase 4 -> aristas inversas
    template <class DistanciaEntre>
    void fase4_agregar_inversas(DistanciaEntre& distancia_entre)
    {
        std::vector<std::vector<slot_t>> con_inversas = vecinos_;
        for (slot_t p = 0; p != cantidad_slots_; ++p)
        {
            for (slot_t u : vecinos_[p]) {
                auto& lista_u = con_inversas[u];
                if (std::find(lista_u.begin(), lista_u.end(), p) == lista_u.end()) {
                    lista_u.push_back(p);
                    ++estadisticas_.aristas_inversas_agregadas;
                }
            }
        }

        for (slot_t u = 0; u != cantidad_slots_; ++u) {
            if (con_inversas[u].size() <= parametros_.max_vecinos) continue;

            std::vector<Neighbor> candidatos;
            for (slot_t vecino : con_inversas[u])
                candidatos.push_back( { distancia_entre(u, vecino), vecino }) ;
            
            std::sort(candidatos.begin(), candidatos.end(), closer);
            con_inversas[u].clear();
            
            for (Neighbor const& elegido : poda_adaptativa(candidatos, parametros_, distancia_entre))
                con_inversas[u].push_back(elegido.slot);
            ++estadisticas_.nodos_repodados;
        }
        vecinos_ = std::move(con_inversas);
    }

    // fase 5 -> conectar islas
    template <class Distancias>
    void fase5_asegurar_conectividad(Distancias const& distancias) {
        std::vector<char> alcanzado(cantidad_slots_, 0);
        marcar_alcanzables(entrada_, alcanzado);
        for (slot_t huerfano = 0; huerfano != cantidad_slots_; ++huerfano) {
            if (alcanzado[huerfano])
                continue;
            auto const distancia_al_huerfano = [&](slot_t slot) { return distancias.entre(huerfano, slot); };
            std::vector<Neighbor> const visitados = recorrer_beam(
                entrada_, parametros_.ancho_beam_construccion, vecinos_, distancia_al_huerfano);

            Neighbor const padre = *std::min_element(visitados.begin(), visitados.end(), closer);

            vecinos_[padre.slot].push_back(huerfano);

            marcar_alcanzables(huerfano, alcanzado);
            ++estadisticas_.nodos_reconectados;
        }
    }

    void marcar_alcanzables(slot_t inicio, std::vector<char>& alcanzado) const {
        std::vector<slot_t> por_visitar{inicio};
        alcanzado[inicio] = 1;
        while (!por_visitar.empty()) {
            slot_t const actual = por_visitar.back();
            por_visitar.pop_back();
            for (slot_t vecino : vecinos_[actual])
                if (!alcanzado[vecino]) {
                    alcanzado[vecino] = 1;
                    por_visitar.push_back(vecino);
                }
        }
    }

    // lista del nodo v; vacia si v se agrego despues del ultimo build()
    std::vector<slot_t> const& lista_de(slot_t v) const {
        static std::vector<slot_t> const vacia;
        return v < vecinos_.size() ? vecinos_[v] : vacia;
    }

    std::size_t cantidad_aristas() const noexcept {
        std::size_t total = 0;
        for (auto const& lista : vecinos_)
            total += lista.size();
        return total;
    }

    double grado_promedio() const {
        std::size_t total = 0;
        for (auto const& lista : vecinos_)
            total += lista.size();
        return cantidad_slots_ ? (double)(total) / (double)(cantidad_slots_) : 0.0;
    }

    ParametrosAlphaCG parametros_;
    slot_t cantidad_slots_ = 0;                     // nodos registrados (vivos o no)
    std::vector< std::vector<slot_t> > vecinos_;    // grafo final: slot -> vecinos de salida
    slot_t entrada_ = 0;                            // nodo donde empieza toda busqueda
    bool grafo_al_dia_ = true;                      // false tras add/update hasta el proximo build()
    Estadisticas estadisticas_;
};

} // namespace eda
