#ifndef ARGUMENTS_HPP
#define ARGUMENTS_HPP

#include <args.hxx>

struct SeedingArguments {
    SeedingArguments(args::ArgumentParser& parser)
        : parser(parser)
        //n{parser, "INT", "Number of strobes [2]", {'n'}}
        , r{parser, "INT",
            "Mean read length. This parameter is estimated from the first 500 "
            "records in each read file. No need to set this explicitly unless you have a reason.", {'r'}}
        , m{parser, "INT",
            "Maximum seed length. Defaults to r - 50. For reasonable values on -l and -u, "
            "the seed length distribution is usually determined by parameters l and u. "
            "Then, this parameter is only active in regions where syncmers are very sparse.", {'m'}}
        , k{parser, "INT", "Strobe length, has to be below 32. [20]", {'k'}}
        , l{parser, "INT", "Lower syncmer offset from k/(k-s+1). Start sample second syncmer k/(k-s+1) + l syncmers downstream [0]", {'l'}}
        , u{parser, "INT", "Upper syncmer offset from k/(k-s+1). End sample second syncmer k/(k-s+1) + u syncmers downstream [7]", {'u'}}
        , c{parser, "INT", "Bitcount length between 2 and 63. [8]", {'c'}}
        , s{parser, "INT",
            "Submer size used for creating syncmers [k-4]. Only even numbers on k-s allowed. "
            "A value of s=k-4 roughly represents w=10 as minimizer window [k-4]. "
            "It is recommended not to change this parameter unless you have a good "
            "understanding of syncmers as it will drastically change the memory usage and "
            "results with non default values.", {'s'}}
    {
    }
    args::ArgumentParser& parser;
    args::ValueFlag<int> r, m, k, l, u, c, s;
};

#endif
