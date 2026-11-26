#pragma once
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <type_traits>
#include <vector>
#include <usearch/index.hpp>
#include <usearch/index_plugins.hpp> 
#include <usearch/alpha_cg/alpha_cg_backend.hpp>

namespace unum {
namespace usearch {

template <typename distance_at, typename key_at, typename compressed_slot_at,
          typename dynamic_allocator_at, typename tape_allocator_at>
class index_alpha_cg_gt {
  public:
    using distance_t = distance_at;
    using vector_key_t = key_at;
    using compressed_slot_t = compressed_slot_at;
    using member_cref_t = member_cref_gt<vector_key_t>;
    using member_ref_t = member_ref_gt<vector_key_t>;
    static_assert(sizeof(compressed_slot_t) == sizeof(eda::slot_t), "slot de USearch != slot de eda");


    struct add_result_t {
        error_t error{};
        std::size_t new_size{};
        std::size_t visited_members{};
        std::size_t computed_distances{};
        std::size_t computed_distances_in_refines{};
        std::size_t computed_distances_in_reverse_refines{};
        compressed_slot_t slot{};

        explicit operator bool() const noexcept { return !error; }
        add_result_t failed(error_t message) noexcept {
            error = std::move(message);
            return std::move(*this);
        }
    };

    /// @brief  Describes a matched search result, augmenting `member_cref_t`
    ///         contents with `distance` to the query object.
    struct match_t {
        member_cref_t member;
        distance_t distance;

        inline match_t() noexcept : member({nullptr, 0}), distance((std::numeric_limits<distance_t>::max)()) {}

        inline match_t(member_cref_t member, distance_t distance) noexcept : member(member), distance(distance) {}

        inline match_t(match_t&& other) noexcept
            : member({other.member.key.ptr(), other.member.slot}), distance(other.distance) {}

        inline match_t(match_t const& other) noexcept
            : member({other.member.key.ptr(), other.member.slot}), distance(other.distance) {}

        inline match_t& operator=(match_t const& other) noexcept {
            member.key.reset(other.member.key.ptr());
            member.slot = other.member.slot;
            distance = other.distance;
            return *this;
        }

        inline match_t& operator=(match_t&& other) noexcept {
            member.key.reset(other.member.key.ptr());
            member.slot = other.member.slot;
            distance = other.distance;
            return *this;
        }
    };

    struct cluster_result_t {
        error_t error{};
        std::size_t visited_members{};
        std::size_t computed_distances{};
        match_t cluster{};

        explicit operator bool() const noexcept { return !error; }
        cluster_result_t failed(error_t message) noexcept {
            error = std::move(message);
            return std::move(*this);
        }
    };

    struct stats_t {
        std::size_t nodes{};
        std::size_t edges{};
        std::size_t max_edges{};
        std::size_t allocated_bytes{};
    };

    using neighbors_view_t = span_gt<compressed_slot_t const>;

    template <typename ref_at, typename index_at> class member_iterator_gt {
        using ref_t = ref_at;
        using index_t = index_at;

        friend class index_alpha_cg_gt;
        member_iterator_gt() noexcept {}
        member_iterator_gt(index_t* index, compressed_slot_t slot) noexcept : index_(index), slot_(slot) {}

        index_t* index_{};
        compressed_slot_t slot_{};

      public:
        using iterator_category = std::random_access_iterator_tag;
        using value_type = ref_t;
        using difference_type = std::ptrdiff_t;
        using pointer = void;
        using reference = ref_t;

        reference operator*() const noexcept { return index_->at(slot_); }
        vector_key_t key() const noexcept { return index_->at(slot_).key; }

        friend inline compressed_slot_t get_slot(member_iterator_gt const& it) noexcept { return it.slot_; }
        friend inline vector_key_t get_key(member_iterator_gt const& it) noexcept { return it.key(); }

        // clang-format off
        member_iterator_gt operator++(int) noexcept { member_iterator_gt old(index_, slot_); ++(*this); return old; }
        member_iterator_gt operator--(int) noexcept { member_iterator_gt old(index_, slot_); --(*this); return old; }
        member_iterator_gt operator+(difference_type d) noexcept { return member_iterator_gt(index_, static_cast<compressed_slot_t>(static_cast<std::size_t>(slot_) + d)); }
        member_iterator_gt operator-(difference_type d) noexcept { return member_iterator_gt(index_, static_cast<compressed_slot_t>(static_cast<std::size_t>(slot_) - d)); }
        member_iterator_gt& operator++() noexcept { slot_ = static_cast<compressed_slot_t>(static_cast<std::size_t>(slot_) + 1); return *this; }
        member_iterator_gt& operator--() noexcept { slot_ = static_cast<compressed_slot_t>(static_cast<std::size_t>(slot_) - 1); return *this; }
        member_iterator_gt& operator+=(difference_type d) noexcept { slot_ = static_cast<compressed_slot_t>(static_cast<std::size_t>(slot_) + d); return *this; }
        member_iterator_gt& operator-=(difference_type d) noexcept { slot_ = static_cast<compressed_slot_t>(static_cast<std::size_t>(slot_) - d); return *this; }
        bool operator==(member_iterator_gt const& other) const noexcept { return index_ == other.index_ && slot_ == other.slot_; }
        bool operator!=(member_iterator_gt const& other) const noexcept { return index_ != other.index_ || slot_ != other.slot_; }
        // clang-format on
    };

    using member_iterator_t = member_iterator_gt<member_ref_t, index_alpha_cg_gt>;
    using member_citerator_t = member_iterator_gt<member_cref_t, index_alpha_cg_gt const>;


    class search_result_t {
        std::vector<eda::Neighbor> encontrados_; 
        vector_key_t const* claves_ = nullptr;
        friend class index_alpha_cg_gt;

      public:
        std::size_t count{};              
        std::size_t visited_members{};   
        std::size_t computed_distances{}; 
        error_t error{};

        search_result_t() noexcept {}
        search_result_t(search_result_t&&) = default;
        search_result_t& operator=(search_result_t&&) = default;

        explicit operator bool() const noexcept { return !error; }
        search_result_t failed(error_t message) noexcept {
            error = std::move(message);
            return std::move(*this);
        }

        operator std::size_t() const noexcept { return count; }
        std::size_t size() const noexcept { return count; }
        bool empty() const noexcept { return !count; }

        match_t at(std::size_t i) const noexcept {
            eda::Neighbor const& vecino = encontrados_[i];
            byte_t* clave = const_cast<byte_t*>(reinterpret_cast<byte_t const*>(claves_ + vecino.slot));
            return {member_cref_t{clave, vecino.slot}, vecino.distancia};
        }
        match_t operator[](std::size_t i) const noexcept { return at(i); }
        match_t front() const noexcept { return at(0); }
        match_t back() const noexcept { return at(count - 1); }

        bool contains(vector_key_t key) const noexcept {
            for (std::size_t i = 0; i != count; ++i)
                if (claves_[encontrados_[i].slot] == key)
                    return true;
            return false;
        }

        std::size_t dump_to(vector_key_t* keys, distance_t* distances) const noexcept {
            for (std::size_t i = 0; i != count; ++i) {
                keys[i] = claves_[encontrados_[i].slot];
                distances[i] = encontrados_[i].distancia;
            }
            return count;
        }
        std::size_t dump_to(vector_key_t* keys) const noexcept {
            for (std::size_t i = 0; i != count; ++i)
                keys[i] = claves_[encontrados_[i].slot];
            return count;
        }
        std::size_t dump_to(vector_key_t* keys, distance_t* distances, std::size_t capacity) const noexcept {
            std::size_t const llenos = (std::min)(count, capacity);
            dump_to(keys, distances);
            for (std::size_t i = llenos; i != capacity; ++i) {
                keys[i] = vector_key_t{};
                distances[i] = (std::numeric_limits<distance_t>::max)();
            }
            return llenos;
        }
    };


    explicit index_alpha_cg_gt(index_config_t config = {}) : grafo_(parametros_desde(config)) {}

    std::size_t size() const noexcept { return grafo_.size(); }
    std::size_t capacity() const noexcept { return limits_.members; }
    std::size_t connectivity() const noexcept { return grafo_.parametros().max_vecinos; }
    std::size_t max_level() const noexcept { return 0; } // una sola capa
    index_limits_t const& limits() const noexcept { return limits_; }
    bool is_immutable() const noexcept { return false; }

    bool try_reserve(index_limits_t limits) noexcept {
        limits.members = (std::max)(limits.members, limits_.members); // nunca achicar
        try {
            claves_.resize(limits.members);
            grafo_.reserve(limits.members);
        } catch (...) {
            return false;
        }
        limits_ = limits;
        return true;
    }
    bool reserve(index_limits_t limits) noexcept { return try_reserve(limits); }

    void clear() noexcept { grafo_.clear(); }
    void reset() noexcept {
        grafo_.clear();
        claves_.clear();
        claves_.shrink_to_fit();
        limits_ = index_limits_t{unreserved};
    }

    member_ref_t at(compressed_slot_t slot) noexcept {
        return {misaligned_ref_gt<vector_key_t>(reinterpret_cast<byte_t*>(&claves_[slot])), slot};
    }
    member_cref_t at(compressed_slot_t slot) const noexcept {
        byte_t* clave = const_cast<byte_t*>(reinterpret_cast<byte_t const*>(&claves_[slot]));
        return {misaligned_ref_gt<vector_key_t const>(clave), slot};
    }
    member_iterator_t iterator_at(compressed_slot_t slot) noexcept { return {this, slot}; }


    template <typename value_at, typename metric_at, typename callback_at = dummy_callback_t,
              typename prefetch_at = dummy_prefetch_t>
    add_result_t add(vector_key_t key, value_at&&, metric_at&&, index_update_config_t = {},
                     callback_at&& callback = callback_at{}, prefetch_at&& = prefetch_at{}) {
        add_result_t result;
        if (size() >= capacity())
            return result.failed("Reserve capacity ahead of insertions!"); // si no, USearch escribe fuera
        compressed_slot_t const slot = static_cast<compressed_slot_t>(size());
        claves_[slot] = key;
        grafo_.add(static_cast<eda::slot_t>(slot));
        result.new_size = size();
        result.slot = slot;
        callback(at(slot)); 
        return result;
    }

    template <typename value_at, typename metric_at, typename callback_at = dummy_callback_t,
              typename prefetch_at = dummy_prefetch_t>
    add_result_t update(member_iterator_t iterator, vector_key_t key, value_at&&, metric_at&&,
                        index_update_config_t = {}, callback_at&& callback = callback_at{},
                        prefetch_at&& = prefetch_at{}) {
        add_result_t result;
        compressed_slot_t const slot = get_slot(iterator);
        claves_[slot] = key;
        grafo_.update(static_cast<eda::slot_t>(slot));
        result.new_size = size();
        result.slot = slot;
        callback(at(slot));
        return result;
    }


    template <typename metric_at>
    error_t build(metric_at&& metric, metric_kind_t tipo, std::size_t /*expansion*/) {
        if (tipo != metric_kind_t::l2sq_k && tipo != metric_kind_t::cos_k)
            return "alpha-CG necesita la metrica l2sq o cos (la poda usa la desigualdad triangular)";
        distancias_t<typename std::decay<metric_at>::type> const distancias{&metric, this, tipo};
        grafo_.build(distancias);
        return {};
    }


    template <typename value_at, typename metric_at, typename predicate_at = dummy_predicate_t,
              typename prefetch_at = dummy_prefetch_t>
    search_result_t search(value_at&& query, std::size_t wanted, metric_at&& metric,
                           index_search_config_t config = {}, predicate_at&& predicate = predicate_at{},
                           prefetch_at&& = prefetch_at{}) const {
        search_result_t result;
        result.claves_ = claves_.data();
        if (!wanted || !size())
            return result;

        distancias_t<typename std::decay<metric_at>::type> const distancias{&metric, this, metric_kind_t::unknown_k};
        auto permitir = [&](eda::slot_t slot) { return predicate(at(slot)); }; // filtra borrados y filtro del usuario
        eda::ContadoresBusqueda contadores;

        if (config.exact) {
            for (eda::slot_t slot = 0; slot != size(); ++slot)
                if (permitir(slot))
                    result.encontrados_.push_back({distancias.desde(query, slot), slot});
            contadores.distancias_calculadas = result.encontrados_.size();
            std::size_t const cuantos = (std::min)(wanted, result.encontrados_.size());
            std::partial_sort(result.encontrados_.begin(), result.encontrados_.begin() + cuantos,
                              result.encontrados_.end(), eda::closer);
            result.encontrados_.resize(cuantos);
        } else {
            result.encontrados_ = grafo_.search(query, wanted, config.expansion, distancias, permitir, &contadores);
        }

        result.count = result.encontrados_.size();
        result.computed_distances = contadores.distancias_calculadas;
        result.visited_members = contadores.nodos_expandidos;
        return result;
    }

 
    std::size_t serialized_length() const noexcept {
        return sizeof(std::uint64_t) + size() * sizeof(vector_key_t) + grafo_.serialized_length();
    }

    template <typename output_callback_at, typename progress_at = dummy_progress_t>
    serialization_result_t save_to_stream(output_callback_at&& output, progress_at&& progress = {}) const {
        serialization_result_t result;
        std::uint64_t const n = size();
        if (!output(&n, sizeof(n)))
            return result.failed("No se pudo escribir la cantidad de claves");
        if (n && !output(claves_.data(), static_cast<std::size_t>(n) * sizeof(vector_key_t)))
            return result.failed("No se pudieron escribir las claves");
        eda::ResultadoSerializacion const grafo = grafo_.save_to_stream(output, progress);
        if (!grafo)
            return result.failed(grafo.error);
        return result;
    }

    template <typename input_callback_at, typename progress_at = dummy_progress_t>
    serialization_result_t load_from_stream(input_callback_at&& input, progress_at&& progress = {}) {
        serialization_result_t result;
        reset();
        std::uint64_t n = 0;
        if (!input(&n, sizeof(n)))
            return result.failed("No se pudo leer la cantidad de claves");
        std::vector<vector_key_t> claves(static_cast<std::size_t>(n));
        if (n && !input(claves.data(), static_cast<std::size_t>(n) * sizeof(vector_key_t)))
            return result.failed("No se pudieron leer las claves");
        eda::ResultadoSerializacion const grafo = grafo_.load_from_stream(input, progress);
        if (!grafo)
            return result.failed(grafo.error);
        if (grafo_.size() != n) {
            reset();
            return result.failed("Archivo corrupto: claves y grafo no coinciden");
        }
        claves_ = std::move(claves);
        limits_.members = static_cast<std::size_t>(n);
        return result;
    }

    template <typename progress_at = dummy_progress_t>
    serialization_result_t view(memory_mapped_file_t, std::size_t = 0, progress_at&& = {}) {
        return serialization_result_t{}.failed("alpha-CG: view() todavia no esta soportado, usa load()");
    }

  private:
    template <typename metric_at> struct distancias_t {
        using consulta_t = byte_t const*;
        metric_at const* metrica;
        index_alpha_cg_gt const* indice;
        metric_kind_t tipo;

        eda::distance_t entre(eda::slot_t a, eda::slot_t b) const { return (*metrica)(indice->at(a), indice->at(b)); }
        eda::distance_t desde(consulta_t q, eda::slot_t s) const { return (*metrica)(q, indice->at(s)); }
        eda::distance_t real(eda::distance_t d) const {
            return tipo == metric_kind_t::cos_k ? std::sqrt(2 * (d > 0 ? d : 0)) // coseno -> cuerda
                                                : std::sqrt(d);                  // l2sq -> l2
        }
    };

    static eda::ParametrosAlphaCG parametros_desde(index_config_t const& config) {
        eda::ParametrosAlphaCG parametros;
        parametros.max_vecinos = config.connectivity;
        return parametros;
    }

    eda::AlphaCGBackend grafo_;
    std::vector<vector_key_t> claves_;
    index_limits_t limits_{unreserved};
};

} // namespace usearch
} // namespace unum
