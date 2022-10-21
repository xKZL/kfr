/**
 * KFR (http://kfrlib.com)
 * Copyright (C) 2016  D Levin
 * See LICENSE.txt for details
 */

#include <kfr/base.hpp>
#include <kfr/dsp.hpp>
#include <kfr/io.hpp>

using namespace kfr;

std::string to_string2(double x)
{
    std::ostringstream ss;
    ss << x;
    return ss.str();
}

int main()
{
    println(library_version());

    int rp = 1;
    int rs = 6;
    constexpr size_t maxorder = 32;
    const char* k[maxorder] = // k coefficients computed with rp = 1 & rs = 6 
    {                          // (k depends on z and p.. so this should be enought in a first time to check the validity)
        "0.89125093813374556273032567332848", // k = 0
        "1.96522672836027156861860021308530", // k = 1
        "0.50118723362727279901918109317194", // k = 2
        "0.98794734469230960360874860270997", // k = 3
        "0.50118723362727934933502638159553", // k = 4
        "0.97538440598394993141795339397504", // k = 5
        "0.50118723362745221105996051846887", // k = 6
        "0.97511086984454220516482791936141", // k = 7
        "0.50118723362753969663430098080426", // k = 8
        "0.97510485666884338940008092322387", // k = 9
        "0.50118723375847462619958605500869", // k = 10
        "0.97510467362126196366745034538326", // k = 11
        "0.50118498126436306083775207298459", // k = 12
        "0.97482310868854160634811023555812", // k = 13
        "0.49486370015702924041178789593687", // k = 14
        "0.83465109606081722137815859241528", // k = 15
        "0.30887740092250570711485124775209", // k = 16 
        "0.48361023745374215332404332912120", // k = 17 
        "0.14515935341003305403262402251130", // k = 18 
        "0.26007573487724366945172960186028", // k = 19 
        "0.06214721652055405637371521265777", // k = 20
        "0.12856256785699293754277050538803", // k = 21
        "0.02471231445727731929062898075244", // k = 22
        "0.05842982129840767341333318540819", // k = 23
        "0.00919341192048158013794267873209", // k = 24
        "0.02453838324842558954452798047896", // k = 25
        "0.00321322302422841340682757582670", // k = 26
        "0.00957488847253192346120620470629", // k = 27
        "0.00105876470543555894175680176517", // k = 28
        "0.00348895553891140223004563303277", // k = 29
        "0.00032993230118648272106499086398", // k = 30
        "0.00119263844932737412413148447854", // k = 31
    };

    for(int i = 0; i < 32; i++) 
    {
        zpk<fbase> filt = kfr::elliptic<fbase>(i, rp, rs);
            
        char buffer[35];
        snprintf(buffer, sizeof(buffer), "%.33f", filt.k);

        std::cout << "Test #" << i << (i < 10 ? " " : "") << ": "<< k[i] << " (python SciPy/stdc++)";
        printf(" <> %s = ", buffer);

        int epsilon = 0;
        bool result = true;
        for(int n = 0; n < 35; n++) {

            result &= (k[i][n] == buffer[n]);
            if(!result) {
                epsilon = n;
                break;
            }
        }

        std::cout << result << " (1e-"<<epsilon<<")";

        std::cout << std::endl;
    }

    return 0;
}
