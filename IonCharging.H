/* SPDX-License-Identifier: GPL-3.0-or-later */
#ifndef IonCharging_H
#define IonCharging_H

#include "CloudFunctionObject.H"
#include "OFstream.H"
#include "OSspecific.H"

#include <cmath>
#include <limits>
#include <map>
#include <memory>

extern "C"
{
#include "charging.h"
}

namespace Foam
{

template<class Cloud>
class IonCharging : public CloudFunctionObject<Cloud>
{
    using Parcel = typename Cloud::parcelType;
    using TrackingData = typename Parcel::trackingData;
    using Base = CloudFunctionObject<Cloud>;

    struct FreeTable
    {
        void operator()(EtaTable* t) const
        {
            eta_table_free(t);
            delete t;
        }
    };

    using Table = std::unique_ptr<EtaTable, FreeTable>;
    using Tables = std::map<scalar, Table>;

    const volScalarField& nP_;
    const volScalarField& nN_;
    const label qmin_;
    const label qmax_;
    std::shared_ptr<Tables> tables_;

    label checkedCharge(const scalar charge) const
    {
        if
        (
            !std::isfinite(charge) || charge != std::trunc(charge)
         || charge < qmin_ || charge > qmax_
        )
        {
            FatalErrorInFunction
                << "Parcel charge " << charge << " must be an integer in ["
                << qmin_ << ", " << qmax_ << "]." << exit(FatalError);
        }
        return label(charge);
    }

    const EtaTable& tableFor(const scalar diameter) const
    {
        const scalar tolerance = 1e-9*diameter;
        const auto iter = tables_->lower_bound(diameter - tolerance);
        if (iter == tables_->end() || mag(iter->first - diameter) > tolerance)
        {
            FatalErrorInFunction
                << "No precomputed ion collision table for diameter "
                << diameter << " m. Use fixedValue injection diameters."
                << exit(FatalError);
        }
        return *iter->second;
    }

    void writeTables() const
    {
        if (!Pstream::master())
        {
            return;
        }
        const fileName output
        (
            this->owner().mesh().time().globalPath()
           /"postProcessing"/"ionCharging"/this->owner().name()
        );
        mkDir(output);
        OFstream csv(output/"Collision_Coefficients.csv");
        csv.precision(16);
        csv << "diameter_m,q,betaP_m3_s,betaN_m3_s,positive_model,negative_model,"
            << "positive_validity,negative_validity" << nl;
        for (const auto& item : *tables_)
        {
            const EtaTable& t = *item.second;
            for (int i = 0; i < t.m; ++i)
            {
                csv << item.first << ',' << i + t.qmin
                    << ',' << t.eta_plus[i] << ',' << t.eta_minus[i]
                    << ',' << collision_model_name(t.diag_plus[i].model_used)
                    << ',' << collision_model_name(t.diag_minus[i].model_used)
                    << ',' << collision_validity_name(&t.diag_plus[i])
                    << ',' << collision_validity_name(&t.diag_minus[i]) << nl;
            }
        }
        Info<< "Wrote ion collision coefficients to " << csv.name() << nl;
    }

public:

    IonCharging
    (
        Cloud& cloud,
        const volScalarField& nP,
        const volScalarField& nN,
        const dictionary& gas
    )
    :
        Base(cloud),
        nP_(nP),
        nN_(nN),
        qmin_(cloud.particleProperties().subDict("constantProperties")
            .template get<label>("qmin")),
        qmax_(cloud.particleProperties().subDict("constantProperties")
            .template get<label>("qmax")),
        tables_(std::make_shared<Tables>())
    {
        const dictionary& props =
            cloud.particleProperties().subDict("constantProperties");
        const scalar eps = props.get<scalar>("epsilonR");
        const dimensionSet mobilityDimensions(-1, 0, 2, 0, 0, 1, 0);
        const scalar muP = dimensionedScalar("muP", mobilityDimensions, gas).value();
        const scalar muN = dimensionedScalar("muN", mobilityDimensions, gas).value();
        const scalar temperature = dimensionedScalar("T", dimTemperature, gas).value();

        if
        (
            qmax_ < qmin_ || int64_t(qmax_) - int64_t(qmin_) >= std::numeric_limits<int>::max()
         || !std::isfinite(eps) || eps < 1
         || !std::isfinite(muP) || muP <= 0
         || !std::isfinite(muN) || muN <= 0
         || !std::isfinite(temperature) || temperature <= 0
        )
        {
            FatalIOErrorInFunction(props)
                << "Invalid charge bounds, relative permittivity, gas temperature, "
                << "or ion mobility." << exit(FatalIOError);
        }

        const dictionary& subModels = cloud.particleProperties().subDict("subModels");
        const dictionary& injections = subModels.subDict("injectionModels");
        for (const entry& model : injections)
        {
            const dictionary& injection = model.dict();
            const word modelType(injection.get<word>("type"));
            if (modelType != "manualInjection" && modelType != "patchInjection")
            {
                FatalIOErrorInFunction(injection)
                    << "Ion charging needs manualInjection or patchInjection "
                    << "with a fixedValue sizeDistribution." << exit(FatalIOError);
            }
            checkedCharge(injection.getOrDefault<scalar>("parcelCharge", 0));
            const dictionary& sizes = injection.subDict("sizeDistribution");
            if (sizes.get<word>("type") != "fixedValue")
            {
                FatalIOErrorInFunction(sizes)
                    << "Use one fixedValue injection entry per particle diameter "
                    << "to precompute all collision coefficients." << exit(FatalIOError);
            }
            const scalar d = sizes.subDict("fixedValueDistribution").get<scalar>("value");
            if (!std::isfinite(d) || d <= 0)
            {
                FatalIOErrorInFunction(sizes)
                    << "Particle diameter must be positive." << exit(FatalIOError);
            }
            if (tables_->count(d))
            {
                continue;
            }
            Table table(new EtaTable{});
            if (eta_table_build(table.get(), d, temperature, qmin_, qmax_, eps, muN, muP))
            {
                FatalErrorInFunction
                    << "Cannot build ion collision coefficients for diameter " << d
                    << " m." << exit(FatalError);
            }
            tables_->emplace(d, std::move(table));
        }
        if (tables_->empty())
        {
            FatalIOErrorInFunction(injections)
                << "At least one injection entry is required." << exit(FatalIOError);
        }
        for (const auto& parcel : cloud)
        {
            checkedCharge(parcel.chargeNumber());
            tableFor(parcel.d());
        }
        writeTables();
    }

    autoPtr<Base> clone() const override
    {
        return autoPtr<Base>(new IonCharging(*this));
    }

    void postEvolve(const TrackingData&) override {}

    bool postMove
    (
        Parcel& parcel,
        const scalar dt,
        const point&,
        const TrackingData& td
    ) override
    {
        if (parcel.active() && dt > 0)
        {
            label charge = checkedCharge(parcel.chargeNumber());
            const EtaTable& table = tableFor(parcel.d());
            const scalar positive = max(scalar(0), nP_[parcel.cell()]);
            const scalar negative = max(scalar(0), nN_[parcel.cell()]);
            auto& rng = this->owner().rndGen();
            scalar remaining = dt;
            while (remaining > 0)
            {
                const label i = charge - qmin_;
                const scalar plus = charge < qmax_ ? table.eta_plus[i]*positive : 0;
                const scalar minus = charge > qmin_ ? table.eta_minus[i]*negative : 0;
                const scalar total = plus + minus;
                if (!std::isfinite(total))
                {
                    FatalErrorInFunction
                        << "Non-finite ion attachment rate." << exit(FatalError);
                }
                if (total <= 0)
                {
                    break;
                }
                const scalar hazard = -std::log
                (
                    max(VSMALL, rng.template sample01<scalar>())
                );
                // Check for an event before dividing by a potentially tiny rate.
                if (hazard >= total*remaining)
                {
                    break;
                }
                const scalar waiting = hazard/total;
                const scalar next = remaining - waiting;
                if (next >= remaining)
                {
                    FatalErrorInFunction
                        << "Ion waiting time is below numerical resolution."
                        << exit(FatalError);
                }
                remaining = next;
                charge += (rng.template sample01<scalar>()*total < minus) ? -1 : 1;
            }
            parcel.chargeNumber() = charge;
        }
        return td.keepParticle;
    }
};

}
#endif
