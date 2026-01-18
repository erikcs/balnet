#pragma once
#include <adelie_core/glm/glm_base.hpp>

#ifndef GLM_CBPS_TP
#define GLM_CBPS_TP \
    template <class ValueType>
#endif
#ifndef GLM_CBPS
#define GLM_CBPS \
    GlmCBPS<ValueType>
#endif

namespace adelie_core {
namespace glm {

template <class ValueType>
class GlmCBPS: public GlmBase<ValueType>
{
public:
    using base_t = GlmBase<ValueType>;
    using typename base_t::value_t;
    using typename base_t::vec_value_t;
    using base_t::y;
    using base_t::weights;

    const value_t target_scale;

    explicit GlmCBPS(
        const Eigen::Ref<const vec_value_t>& y,
        const Eigen::Ref<const vec_value_t>& weights,
        value_t target_scale
    );

    ADELIE_CORE_GLM_PURE_OVERRIDE_DECL
};

} // namespace glm
} // namespace adelie_core
