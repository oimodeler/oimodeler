# -*- coding: utf-8 -*-
"""
Created on Fri Oct 18 12:15:03 2024

@author: ame
"""

import numpy as np
from astroquery.vizier import Vizier

import oimodeler as oim

JSDC_CAT = "II/346/jsdc_v2"
BAND_NAMES = np.array(["B", "V", "R", "I", "J", "H", "K", "L", "M", "N"])
WL_BANDS = (
    np.array([445, 551, 658, 806, 1220, 1630, 2190, 3450, 4750, 10500]) * 1e-9
)


def getJSDCDiameter(name, bandOrWl):
    res = Vizier.query_object(name, catalog=JSDC_CAT)

    if type(bandOrWl) == str:
        band = bandOrWl
    else:
        band = BAND_NAMES[np.argmin(np.abs(WL_BANDS - bandOrWl))]

    if len(res) != 0:
        D = res[0][f"UDD{band}"].data[0]
        eD = res[0]["e_LDD"].data[0]

        return D, eD, band
    else:
        return None, None


def _phaser(phase):
    return np.exp(complex(0, 1) * np.deg2rad(phase))


def _phase(phaser):
    return np.rad2deg(np.angle(phaser))


def oimCalibrate(sci, cal, D=None, eD=None):

    data_sci = oim.oimData(sci)
    data_cal = oim.oimData(cal)

    if not (D):
        calname = data_cal.data[0]["OI_TARGET"].data["TARGET"][0]
        wlmean = np.mean(data_cal.data[0]["OI_WAVELENGTH"].data["EFF_WAVE"])
        print(
            f"No diameter given for cal {calname} at wl={wlmean*1e6:.1f}um,"
            " searching in JSDC..."
        )
        D, eD, band = getJSDCDiameter(calname, wlmean)
        print(f" UDD{band} = {D:.2f} ± {eD:.2f} mas")

    ud = oim.oimUD(d=D)
    mud = oim.oimModel(ud)
    sim = oim.oimSimulator(data_cal, mud)
    sim.compute(computeSimulatedData=True)

    print(data_sci.info())
    for iarr in range(len(data_sci.data[0])):
        arrname = data_sci.data[0][iarr].name

        if arrname == "OI_VIS2":

            # -------------------------VIS2DATA---------------------------------
            data_cal.data[0]["OI_VIS2"].data[
                "VIS2DATA"
            ] /= sim.simulatedData.data[0]["OI_VIS2"].data["VIS2DATA"]
            data_sci.data[0]["OI_VIS2"].data["VIS2DATA"] /= data_cal.data[0][
                "OI_VIS2"
            ].data["VIS2DATA"]

            # TODO: Put here the code for error computation on VIS2DATA

        if arrname == "OI_VIS":

            # ---------------------------VISAMP---------------------------------
            if (
                True
            ):  # TODO put here a check on the presence of VISAMP data and the kind of DATA
                data_cal.data[0]["OI_VIS"].data[
                    "VISAMP"
                ] /= sim.simulatedData.data[0]["OI_VIS"].data["VISAMP"]
                data_sci.data[0]["OI_VIS"].data["VISAMP"] /= data_cal.data[0][
                    "OI_VIS"
                ].data["VISAMP"]

            # TODO: Put here the code for error computation on VISAMP

            #
            #
            #

            # --------------------------VISPHI---------------------------------
            phaser_model = _phaser(
                sim.simulatedData.data[0]["OI_VIS"].data["VISPHI"]
            )
            phaser_cal = _phaser(data_cal.data[0]["OI_VIS"].data["VISPHI"])
            phaser_sci = _phaser(data_sci.data[0]["OI_VIS"].data["VISPHI"])

            phaser_cal *= np.conjugate(phaser_model)
            phaser_sci *= np.conjugate(phaser_cal)

            data_sci.data[0]["OI_VIS"].data["VISPHI"] = _phase(phaser_sci)

            # TODO: Put here the code for error computation on VISPHI

            #
            #
            #

        if arrname == "OI_T3":

            # ----------------------------T3AMP---------------------------------
            if (
                False
            ):  # TODO put here a check on the presence of VISAMP data and the kind of DATA
                data_cal.data[0]["OI_T3"].data[
                    "T3AMP"
                ] /= sim.simulatedData.data[0]["OI_T3"].data["T3AMP"]
                data_sci.data[0]["OI_T3"].data["T3AMP"] /= data_cal.data[0][
                    "OI_T3"
                ].data["T3AMP"]

            # ----------------------------T3PHI---------------------------------
            phaser_model = _phaser(
                sim.simulatedData.data[0]["OI_T3"].data["T3PHI"]
            )
            phaser_cal = _phaser(data_cal.data[0]["OI_T3"].data["T3PHI"])
            phaser_sci = _phaser(data_sci.data[0]["OI_T3"].data["T3PHI"])

            phaser_cal *= np.conjugate(phaser_model)
            phaser_sci *= np.conjugate(phaser_cal)

            data_sci.data[0]["OI_T3"].data["T3PHI"] = _phase(phaser_sci)

        if arrname == "OI_FLUX":

            # --------------------------FLUXDATA----------------------------
            data_cal.data[0]["OI_FLUX"].data[
                "FLUXDATA"
            ] /= sim.simulatedData.data[0]["OI_FLUX"].data["FLUXDATA"]
            data_sci.data[0]["OI_FLUX"].data["FLUXDATA"] /= data_cal.data[0][
                "OI_FLUX"
            ].data["FLUXDATA"]

    return data_sci
