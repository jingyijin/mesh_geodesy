#include "interval.hpp"

Interval::Interval(Handle e_, double b0_, double b1_, double d0_, double d1_, double sigma_, bool tau_)
: MxHeapable()
{
    e = e_;
    b0 = b0_; b1 = b1_;
    d0 = d0_; d1 = d1_;
    sigma = sigma_;
    tau = tau_;
    xi = 0.0;
    state = Discovered;
    compute_min_dist();
}

Interval::Interval(Handle e_, double b0_, double b1_, double d0_, double d1_, double sigma_, bool tau_, TraversalState state_)
: MxHeapable()
{
    e = e_;
    b0 = b0_; b1 = b1_;
    d0 = d0_; d1 = d1_;
    sigma = sigma_;
    tau = tau_;
    xi = 0.0;
    state = state_;
    compute_min_dist();
}

Interval::Interval(Handle e_, double b0_, double b1_, double d0_, double d1_, double sigma_, double xi_, bool tau_, TraversalState state_)
: MxHeapable()
{
    e = e_;
    b0 = b0_; b1 = b1_;
    d0 = d0_; d1 = d1_;
    sigma = sigma_;
    tau = tau_;
    xi = xi_;
    state = state_;
    compute_min_dist();
}

double Interval::min_dist() const
{
    return sigma + ((d0 < d1) ? d0 : d1);
}

void Interval::compute_min_dist()
{
    heap_key(-min_dist());
}

bool Interval::does_intersect(Interval* iv) const
{
    return ( FLT(get_b0(), iv->get_b1()) && FGT(get_b1(), iv->get_b0()) );
}

bool Interval::cover_whole(double end0, double end1) const
{
    return ( FEQ(get_b0(), end0) && FEQ(get_b1(), end1) );
}

bool Interval::is_very_small() const
{
    return FEQ(get_b0()-get_b1(), 0.0, 1e-4);
}

Interval* Interval::copy(double new_b0, double new_b1, double new_d0, double new_d1)
{
    return new Interval(e, new_b0, new_b1, new_d0, new_d1, sigma, tau, state);
}

void Interval::set_range(double new_b0, double new_b1, double new_d0, double new_d1)
{
    b0 = new_b0; b1 = new_b1;
    d0 = new_d0; d1 = new_d1;
}

void Interval::merge(Interval* a, Vec2& new_s, double new_sigma, double new_xi)
{
    this->set_b0(a->get_b0());
    this->set_d0(norm(new_s-Vec2(a->get_b0(), 0)));
    this->set_d1(norm(new_s-Vec2(this->get_b1(), 0)));
    this->set_sigma(new_sigma);
    this->set_xi(new_xi);
}