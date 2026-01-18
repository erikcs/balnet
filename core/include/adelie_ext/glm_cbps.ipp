/* Calibration loss for mu_1 = E[Y(1)] with logit link.
 * The loss is:
 * l(eta) = sum_{i=1}^{n} w_i (y * exp(-eta) + (1 - y) * eta), where y is 0 or 1
 * with (negative) gradient: -w_i * (-y * exp(-eta) + (1 - y))
 * and raw hessian: y * exp(-eta)
 * We get the mu_0 loss by using the same construction on inverted
 * outcomes y' = 1 - y. The mu_0 loss is used for the ATT, and we
 * scale the gradient with target_scale = n / n_0 to get printed
 * path metrics on the same scale.
 */
#pragma once
#define _USE_MATH_DEFINES // for C++
#include <cmath>
#include <unsupported/Eigen/SpecialFunctions>
#include <adelie_ext/glm_cbps.hpp>
#include <adelie_core/util/macros.hpp>

namespace adelie_core {
namespace glm {

GLM_CBPS_TP
GLM_CBPS::GlmCBPS(
    const Eigen::Ref<const vec_value_t>& y,
    const Eigen::Ref<const vec_value_t>& weights,
    value_t target_scale
):
    base_t("cbps", y, weights),
    target_scale(target_scale)
{}

GLM_CBPS_TP
void
GLM_CBPS::gradient(
    const Eigen::Ref<const vec_value_t>& eta,
    Eigen::Ref<vec_value_t> grad
)
{
    base_t::check_gradient(eta, grad);
    grad = target_scale * weights * (y * (-eta).exp() - (1 - y));
}

GLM_CBPS_TP
void
GLM_CBPS::hessian(
    const Eigen::Ref<const vec_value_t>& eta,
    const Eigen::Ref<const vec_value_t>& grad,
    Eigen::Ref<vec_value_t> hess
)
{
    base_t::check_hessian(eta, grad, hess);
    hess = grad + target_scale * weights * (1 - y);
}

GLM_CBPS_TP
typename GLM_CBPS::value_t
GLM_CBPS::loss(
    const Eigen::Ref<const vec_value_t>& eta
)
{
    base_t::check_loss(eta);
    return (target_scale * weights * (y * (-eta).exp() + (1 - y) * eta)).sum();
}

GLM_CBPS_TP
typename GLM_CBPS::value_t
GLM_CBPS::loss_full()
{
    // we don't use loss full for anything in cbps
    value_t loss = 0;
    return loss;
}

GLM_CBPS_TP
void
GLM_CBPS::inv_link(
    const Eigen::Ref<const vec_value_t>& eta,
    Eigen::Ref<vec_value_t> out
)
{
    out = 1 / (1 + (-eta).exp());
}

} // namespace glm
} // namespace adelie_core
