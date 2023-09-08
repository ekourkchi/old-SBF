import sys, os
from astropy import units as u
from astropy.coordinates import SkyCoord
import bs4
from urllib.request import urlopen as uReq
from bs4 import BeautifulSoup as soup
import pandas as pd


def xcmd(cmd, verbose=False):

    if verbose:
        print("\n" + cmd)

    tmp = os.popen(cmd)
    output = ""
    for x in tmp:
        output += x
    if "abort" in output:
        failure = True
    else:
        failure = tmp.close()
    if failure:
        print("execution of %s failed" % cmd)
        print("error is as follows", output)
        sys.exit()
    else:
        return output


def mySoup(my_url):

    uClient = uReq(my_url)
    page_html = uClient.read()
    uClient.close()
    page_soup = soup(page_html, "html.parser")

    return page_soup


def to_degree(RA, Dec):
    c = SkyCoord(RA, Dec)
    return c.ra.degree, c.dec.degree


def fix_name(name):
    return "+".join([x.upper() for x in name.strip().rsplit(" ") if x != ""])


def get_ned_info(name_list):

    query = ""
    for i, name in enumerate(name_list):
        if i == 0:
            query += fix_name(name)
        else:
            query = query + "%0D%0A" + fix_name(name)

    # query = "UGC+12517%0D%0ANGC+2557"
    my_url = "https://ned.ipac.caltech.edu/cgi-bin/gmd?uplist={}&delimiter=bar&NO_LINKS=1&nondb=row_count&nondb=user_name_msg&nondb=user_objname&crosid=objname&position=ra%2Cdec&position=gallon%2Cgallat&position=bhextin&attdat=attned&gphotoms=q_value&gphotoms=q_unc&gphotoms=ned_value&gphotoms=ned_unc&diamdat=ned_maj_dia&distance=avg&distance=stddev_samp".format(
        query
    )
    page_soup = mySoup(my_url)

    rows = page_soup.findAll("pre")[0].text.split("\n")

    columns = None
    output = {}

    for i, row in enumerate(rows):

        cols = row.split("|")
        cl = []
        for c in cols:
            cc = c.strip()
            if cc != "":
                if cc == "Burstein":
                    cc = "Extinction"
                cl.append(cc)

        try:
            if cl[0] == "Row":
                columns = cl
        except:
            pass
        try:
            if int(cl[0]) > 0:
                output[cl[1]] = cl

        except:
            pass

    if columns and output:
        df = pd.DataFrame.from_dict(output, orient="index", columns=columns)
        df["ra_dec"] = df.apply(lambda row: to_degree(row.RA, row.Dec), axis=1)
        return df

    return None


if __name__ == "__main__":

    name_list = ["m31", "n2418"]

    df = get_ned_info(name_list)
    print(df.head())
