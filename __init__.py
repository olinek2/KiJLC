import csv
import os
import pcbnew
import re
import wx
from decimal import Decimal, getcontext
from pathlib import Path

ref_ignore = ["TP", "T", "NT", "REF**", "G", "H"]

# original rotation db from:
# https://github.com/matthewlai/JLCKicadTools/blob/master/jlc_kicad_tools/cpl_rotations_db.csv
rotations = {
    "^R_Array_Convex_": 90,
    "^R_Array_Concave_": 90,
    "^SOT-223": 180,
    "^SOT-23": 180,
    "^SOT-89": 180,
    "^TSOT-23": 180,
    "^SOT-353": 180,
    "^SOT-363": 180,
    "^QFN-": 270,
    "^LQFP-": 270,
    "^TQFP-": 270,
    "^SOP-4_": 0,
    "^SOP-(?!18_)": 270,
    "^TSSOP-": 270,
    "^SSOP-": 270,
    "^DFN-": 270,
    "^SOIC-": 270,
    "^SOP-18_": 0,
    "^VSSOP-8_3.0x3.0mm_P0.65mm": 270,
    "^VSSOP-8_": 180,
    "^VSSOP-10_": 270,
    "^VSON-8_": 270,
    "^TSOP-6": 270,
    "^UDFN-10": 270,
    "^USON-10": 270,
    "^TDSON-8-1": 270,
    "^CP_EIA-": 180,
    "^CP_Elec_8x5.4": 180,
    "^CP_Elec_8x10.5": 180,
    "^CP_Elec_6.3x7.7": 180,
    "^CP_Elec_8x6.7": 180,
    "^CP_Elec_8x10": 180,
    "^CP_Elec_10x10": 180,
    "^(.*?_|V)?QFN-(16|20|24|28|40)(-|_|$)": 270,
    "^Bosch_LGA-8_2x2.5mm_P0.65mm_ClockwisePinNumbering": 90,
    "^PowerPAK_SO-8_Single": 270,
    "^HTSSOP-28-1EP_4.4x9.7mm*": 270,
    "^PUIAudio_SMT_0825_S_4_R*": 270,
    "^USB_C_Receptacle_HRO_TYPE-C-31-M-12*": 180,
    "^ESP32-W": 270,
    "^SOIC127P798X216-8N": -90,
    "^SW_DIP_SPSTx01_Slide_Copal_CHS-01B_W7.62mm_P1.27mm": -180,
    "^BatteryHolder_Keystone_1060_1x2032": -180,
    "^SO-14": -90,
    "^HTSSOP-": 270,
    "^USB_C_Receptacle_XKB_U262-16XN-4BVC11": 180,
    "^Relay_DPDT_Omron_G6K-2F-Y": 270,
    "^RP2040-QFN-56": 270,
    "^TO-277": 90,
    "^LED_WS2812B": 180,
    "^SW_SPST_B3": 90,
    "^Transformer_Ethernet_Pulse_HX0068ANL": 270,
    "^JST_GH_SM": 180,
    "^JST_PH_S": 180,
    "^Diodes_PowerDI3333-8": 270,
}

midpoint_correction = {
    "^PinSocket_1x04_P2.54mm_Vertical": (Decimal(0), Decimal(-3.81)),
    "^PinSocket_1x05_P2.54mm_Vertical": (Decimal(0), Decimal(-5.08)),
    "^XCVR_ESP32-WROVER-B": (Decimal(0), Decimal(0.04)),
    "BarrelJack": (Decimal(-6.5), Decimal(0)),
    "^SW_SPST_HRO": (Decimal(0), Decimal(1.65)),
    "USB_C_Receptacle_HRO_TYPE-C-31-M-12": (Decimal(1.8), Decimal(0.65)),
    "SW_Tactile_SPST_Angled_PTS645Vx58-2LFS": (Decimal(2.2), Decimal(-1)),
}

#
# helper functions from https://docs.python.org/3/library/decimal.html
#

def pi():
    """Compute Pi to the current precision.

    >>> print(pi())
    3.141592653589793238462643383

    """
    getcontext().prec += 2  # extra digits for intermediate steps
    three = Decimal(3)      # substitute "three=3.0" for regular floats
    lasts, t, s, n, na, d, da = 0, three, 3, 1, 0, 0, 24
    while s != lasts:
        lasts = s
        n, na = n+na, na+8
        d, da = d+da, da+32
        t = (t * n) / d
        s += t
    getcontext().prec -= 2
    return +s               # unary plus applies the new precision

def cos(x):
    """Return the cosine of x as measured in radians.

    The Taylor series approximation works best for a small value of x.
    For larger values, first compute x = x % (2 * pi).

    >>> print(cos(Decimal('0.5')))
    0.8775825618903727161162815826
    >>> print(cos(0.5))
    0.87758256189
    >>> print(cos(0.5+0j))
    (0.87758256189+0j)

    """
    getcontext().prec += 2
    i, lasts, s, fact, num, sign = 0, 0, 1, 1, 1, 1
    while s != lasts:
        lasts = s
        i += 2
        fact *= i * (i-1)
        num *= x * x
        sign *= -1
        s += num / fact * sign
    getcontext().prec -= 2
    return +s

def sin(x):
    """Return the sine of x as measured in radians.

    The Taylor series approximation works best for a small value of x.
    For larger values, first compute x = x % (2 * pi).

    >>> print(sin(Decimal('0.5')))
    0.4794255386042030002732879352
    >>> print(sin(0.5))
    0.479425538604
    >>> print(sin(0.5+0j))
    (0.479425538604+0j)

    """
    getcontext().prec += 2
    i, lasts, s, fact, num, sign = 1, 0, x, 1, x, 1
    while s != lasts:
        lasts = s
        i += 2
        fact *= i * (i-1)
        num *= x * x
        sign *= -1
        s += num / fact * sign
    getcontext().prec -= 2
    return +s

class JLCSMTPlugin(pcbnew.ActionPlugin):
    def defaults(self):
        self.name = "Generate JLCSMT Placement Files"
        self.category = "Fabrication Outputs"
        self.description = "Generates the CPL placement files as expected by JLCSMT"
        self.show_toolbar_button = True
        self.icon_file_name = os.path.join(
            os.path.dirname(__file__), 'KiJLC_32x32.png')

    def Run(self):
        board = pcbnew.GetBoard()
#        modules = board.GetModules()
        if hasattr(board, 'GetModules'):
            modules = board.GetModules()
        else:
            modules = board.GetFootprints()
        fn = Path(board.GetFileName()).with_suffix("")

        bot = open("{}_cpl_bot.csv".format(fn), "w", newline='')
        top = open("{}_cpl_top.csv".format(fn), "w", newline='')
        botw = csv.writer(bot, delimiter=',', quotechar='"',
                          quoting=csv.QUOTE_ALL)
        topw = csv.writer(top, delimiter=',', quotechar='"',
                          quoting=csv.QUOTE_ALL)

        hdr = ["Designator", "Mid X", "Mid Y", "Layer", "Rotation"]
        botw.writerow(hdr)
        topw.writerow(hdr)

        for mod in modules:
            skip = False
            ref = mod.GetReference()

            for prefix in ref_ignore:
                if ref.startswith(prefix):
                    skip = True

            if skip:
                continue

            pos = mod.GetPosition()
            rot = mod.GetOrientationDegrees()
            desc = mod.GetDescription()
            layer = board.GetLayerName(mod.GetLayer())
            mid_x = Decimal(pos[0]) / Decimal(1000000)
            mid_y = Decimal(pos[1]) / Decimal(-1000000)
            footprint = str(mod.GetFPID().GetLibItemName())

            print(footprint)

            # some library parts have a different origin than the JLC parts, try to correct it
            for exp in midpoint_correction:
                if re.match(exp, footprint):
                    px, py = midpoint_correction[exp]
                    rad = Decimal(rot) * pi() / Decimal(180)

                    qx = cos(rad) * px - sin(rad) * py
                    qy = sin(rad) * px + cos(rad) * py

                    qx = qx.quantize(Decimal('0.001'))
                    qy = qy.quantize(Decimal('0.001'))

                    print(f"previous midpoint for {footprint} x: {mid_x}, y: {mid_y}; new x: {mid_x + qx}, y: {mid_y + qy}")
                    mid_x += qx
                    mid_y += qy
            
            for exp in rotations:
                if re.match(exp, footprint):
                    new_rot =  (rot + rotations[exp]) % 360
                    print(f"rotating {ref} ({footprint}): prev {rot}, new {new_rot}")
                    rot = new_rot

            x = str(mid_x) + "mm"
            y = str(mid_y) + "mm"

            if layer == "F.Cu":
                topw.writerow([ref, x, y, "top", rot])
            elif layer == "B.Cu":
                botw.writerow([ref, x, y, "bottom", rot])

        bot.close()
        top.close()
        wx.MessageBox("Placement files generated.",
                      'Done', wx.OK | wx.ICON_INFORMATION)


JLCSMTPlugin().register()
