from ctypes import byref, cdll, c_int, c_char_p, c_double, c_float, c_bool
from ctypes import POINTER as pt
import ctypes
import os, sys


class LikeNew:
    def __init__(self, root):

        mylib = cdll.LoadLibrary(
            os.path.join(os.path.dirname(__file__), "lib/likenew6.so")
        )
        self.likenew_ = mylib.likenew_

        self.likenew_.argtypes = [
            c_char_p,
            pt(c_int),
            pt(c_int),
            pt(c_float),
            pt(c_float),
            pt(c_float),
            pt(c_float),
            pt(c_float),
            pt(c_float),
            pt(c_float),
            c_char_p,
            pt(c_int),
            pt(c_bool),
            pt(c_bool),
            pt(c_float),
            pt(c_float),
            pt(c_float),
            pt(c_float),
            pt(c_float),
            pt(c_float),
            pt(c_float),
            pt(c_float),
            pt(c_float),
            pt(c_float),
            pt(c_float),
            pt(c_int),
        ]

        self.root = root

    def run(
        self,
        fname="u12517se.lknj",
        icolor=5,
        secperpixel=0.128,
        fwhm=1.4,
        distance=60,
        kscale=1.2,
        delta=1.4,
        snlim=4.5,
        mlim=22.0,
        mname="u12517sej.ptm6",
        yessoft=False,
        abmags=True,
        verbose=False,
    ):

        fname = self.root + fname
        mname = self.root + mname

        fname = fname.encode("ascii")
        mname = mname.encode("ascii")

        n_fname = len(fname)
        n_mname = len(mname)

        beta = c_float()
        cnorm = c_float()
        cmax = c_float()
        alpha = c_float()
        tot_gc = c_float()
        gamma = c_float()
        gnorm = c_float()
        tnorm = c_float()
        tot_gal = c_float()
        galpersec = c_float()
        tysonempersec = c_float()
        status = c_int()

        self.likenew_(
            fname,
            byref(c_int(n_fname)),
            byref(c_int(icolor)),
            byref(c_float(secperpixel)),
            byref(c_float(fwhm)),
            byref(c_float(distance)),
            byref(c_float(kscale)),
            byref(c_float(delta)),
            byref(c_float(snlim)),
            byref(c_float(mlim)),
            mname,
            byref(c_int(n_mname)),
            byref(c_bool(yessoft)),
            byref(c_bool(abmags)),
            byref(beta),
            byref(cnorm),
            byref(cmax),
            byref(alpha),
            byref(tot_gc),
            byref(gamma),
            byref(gnorm),
            byref(tnorm),
            byref(tot_gal),
            byref(galpersec),
            byref(tysonempersec),
            byref(status),
        )

        out_dict = {}

        # output values from fortran subroutine
        out_dict["beta"] = beta.value
        out_dict["cnorm"] = cnorm.value
        out_dict["cmax"] = cmax.value
        out_dict["alpha"] = alpha.value
        out_dict["tot_gc"] = tot_gc.value
        out_dict["gamma"] = gamma.value
        out_dict["gnorm"] = gnorm.value
        out_dict["tnorm"] = tnorm.value
        out_dict["tot_gal"] = tot_gal.value

        out_dict["galpersec"] = galpersec.value
        out_dict["tysonempersec"] = tysonempersec.value

        # input values of the subroutine
        out_dict["delta"] = delta

        out_dict["status"] = status.value

        if verbose:
            print_out_dict(out_dict)

        return out_dict


def print_out_dict(out_dict):

    print("Beta  =\t      %.4e \t #gxy / total sources" % out_dict["beta"])
    print("")
    print("Cnorm =\t      %.4e \t GC normalization" % out_dict["cnorm"])
    print("Cmax  =\t      %.3f \t\t GC peak magnitude" % out_dict["cmax"])
    print("Delta =\t      %.1f \t\t GC distribution width" % out_dict["delta"])
    print("Alpha =\t      %.3f \t\t GC log slope vs log r" % out_dict["alpha"])
    print("Total # GC =  %.1f" % out_dict["tot_gc"])
    print("")
    print("Gamma =\t      %.5f \t\t Gxy log slope vs m" % out_dict["gamma"])
    print(
        'Gnorm =\t      %.5f \t\t Gxy count / 1/" @ %5.2f'
        % (out_dict["gnorm"], out_dict["galpersec"])
    )
    print(
        "Tnorm =\t      %.5f \t\t Tyson Gxy count @ %5.2f"
        % (out_dict["tnorm"], out_dict["tysonempersec"])
    )
    print("Total # gxy =  %.1f" % out_dict["tot_gal"])


if __name__ == "__main__":

    l = LikeNew()
    l.run()
