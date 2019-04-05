//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include "itensor/util/print_macro.h"
#include "itensor/mps/mpo.h"
#include "itensor/mps/localop.h"

namespace itensor {


using std::istream;
using std::ostream;
using std::cout;
using std::endl;
using std::vector;
using std::find;
using std::pair;
using std::make_pair;
using std::string;

void 
nmultMPO(MPO const& Aorig, 
         MPO const& Borig, 
         MPO& res,
         Args args)
    {
    if(!args.defined("Cutoff")) args.add("Cutoff",1E-14);

    if(length(Aorig) != length(Borig)) Error("nmultMPO(MPO): Mismatched MPO length");
    const int N = length(Borig);

    auto A = Aorig;
    A.position(1);

    MPO B;
    if(&Borig == &Aorig)
        {
        B = A;
        }
    else
        {
        B = Borig;
        B.position(1);
        }

    B.prime();

    res=A;
    auto siA = uniqueIndex(A(1),{B(1),A(2)});
    auto siB = uniqueIndex(B(1),{A(1),B(2)});
    res.ref(1) = ITensor(siA,siB,linkIndex(A,1));

    ITensor clust,nfork;
    for(int i = 1; i < N; ++i)
        {
        if(i == 1) 
            { 
            clust = A(i) * B(i); 
            }
        else       
            { 
            clust = nfork * A(i) * B(i); 
            }

        if(i == N-1) break;

        nfork = ITensor(linkIndex(A,i),linkIndex(B,i),linkIndex(res,i));

        denmatDecomp(clust,res.ref(i),nfork,Fromleft,args);

        auto mid = commonIndex(res(i),nfork,"Link");
        mid.dag();
        auto siA = uniqueIndex(A(i+1),{A(i),A(i+2),B(i+1)});
        auto siB = uniqueIndex(B(i+1),{B(i),B(i+2),A(i+1)});
        res.ref(i+1) = ITensor(mid,siA,siB,rightLinkIndex(res,i+1));
        }

    nfork = clust * A(N) * B(N);

    res.svdBond(N-1,nfork,Fromright, args);
    for(auto i : range1(N))
        {
        if(i < N)
            {
            auto l = linkIndex(res,i);
            res.ref(i).noPrime(l);
            res.ref(i+1).noPrime(l);
            }
        res.ref(i).replaceTags("2","1");
        }
    res.orthogonalize();
    }

//TODO: complete this version that is independent of tag convention
//void 
//nmultMPO(MPO const& Aorig, 
//         MPO const& Borig, 
//         MPO& C,
//         Args args)
//    {
//    if(!args.defined("Cutoff")) args.add("Cutoff",1E-14);
//
//    if(length(Aorig) != length(Borig)) Error("nmultMPO(MPO): Mismatched MPO length");
//    const int N = length(Borig);
//
//    auto A = Aorig;
//    A.position(1);
//
//    MPO B;
//    if(&Borig == &Aorig)
//        {
//        B = A;
//        }
//    else
//        {
//        B = Borig;
//        B.position(1);
//        }
//
//    C = A;
//
//    auto siA = uniqueIndex(A(1),{B(1),A(2)});
//    auto siB = uniqueIndex(B(1),{A(1),B(2)});
//    auto liA = linkIndex(A,1);
//    auto liB = linkIndex(B,1);
//    auto center = A(1) * B(1); 
//
//    auto tagsC = tags(liA);
//
//    ITensor nfork;
//    Index liC;
//    std::tie(C.ref(1),nfork,liC) = denmatDecomp(center,{siA,siB},{liA,liB},Fromleft,{args,"Tags=",tagsC});
//
//    for( auto i = 2; i < N; i++ )
//        {
//        // TODO: use siteInds(A,B,i);
//        siA = uniqueIndex(A(i),{A(i-1),A(i+1),B(i)});
//        siB = uniqueIndex(B(i),{B(i-1),B(i+1),A(i)});
//        liA = linkIndex(A,i);
//        liB = linkIndex(B,i);
//
//        center = nfork * A(i) * B(i); 
//        tagsC = tags(liA);
//        std::tie(C.ref(i),nfork,liC) = denmatDecomp(center,{siA,siB,liC},{liA,liB},Fromleft,{args,"Tags=",tagsC});
//        }
//    C.ref(N) = nfork * A(N) * B(N);
//    C.orthogonalize();
//    }

MPO
nmultMPO(MPO const& A,
         MPO const& B,
         Args args)
  {
  MPO res;
  nmultMPO(A,B,res,args);
  return res;
  }

//
// Define specific applyMPO methods
//

MPS
densityMatrixApplyMPOImpl(MPO const& K,
                          MPS const& x,
                          Args args = Args::global());

void
fitApplyMPOImpl(MPS const& psi,
                MPO const& K,
                MPS & res,
                Args const& args = Args::global());

MPS
applyMPO(MPO const& K,
         MPS const& x,
         Args const& args)
    {
    if( !x ) Error("Error in applyMPO, MPS is uninitialized.");
    if( !K ) Error("Error in applyMPO, MPO is uninitialized.");

    auto method = args.getString("Method","DensityMatrix");

    MPS res;
    if(method == "DensityMatrix")
        {
        res = densityMatrixApplyMPOImpl(K,x,args);
        }
    else if(method == "Fit")
        {
        // Use the input MPS x to be applied as the
        // default starting state
        // TODO: consider using zipUpApplyMPOImpl as 
        // a way to get a better starting state
        auto sites = siteInds(K,x);
        res = replaceSiteInds(x,sites);
        //res = x;
        fitApplyMPOImpl(x,K,res,args);
        }
    else
        {
        Error("applyMPO currently supports the following methods: 'DensityMatrix', 'Fit'");
        }

    return res;
    }


MPS
applyMPO(MPO const& K,
         MPS const& x,
         MPS const& x0,
         Args const& args)
    {
    if( !x ) Error("Error in applyMPO, MPS is uninitialized.");
    if( !K ) Error("Error in applyMPO, MPO is uninitialized.");
    if( !x0 ) Error("Error in applyMPO, guess MPS is uninitialized.");

    auto method = args.getString("Method","Fit");

    MPS res = x0;
    if(method == "DensityMatrix")
        Error("applyMPO method 'DensityMatrix' does not accept an input MPS");
    else if(method == "Fit")
        fitApplyMPOImpl(x,K,res,args);
    else
        Error("applyMPO currently supports the following methods: 'DensityMatrix', 'Fit'");

    return res;
    }

//
// Implement specific applyMPO methods
//


MPS
densityMatrixApplyMPOImpl(MPO const& K,
                          MPS const& psi,
                          Args args)
    {
    if( args.defined("Maxm") )
      {
      if( args.defined("MaxDim") )
        {
        Global::warnDeprecated("Args Maxm and MaxDim are both defined. Maxm is deprecated in favor of MaxDim, MaxDim will be used.");
        }
      else
        {
        Global::warnDeprecated("Arg Maxm is deprecated in favor of MaxDim.");
        args.add("MaxDim",args.getInt("Maxm"));
        }
      }

    auto cutoff = args.getReal("Cutoff",1E-13);
    auto dargs = Args{"Cutoff",cutoff};
    auto maxdim_set = args.defined("MaxDim");
    if(maxdim_set) dargs.add("MaxDim",args.getInt("MaxDim"));
    auto verbose = args.getBool("Verbose",false);
    auto normalize = args.getBool("Normalize",false);

    auto N = length(psi);

    for( auto n : range1(N) )
      {
      if( commonIndex(psi(n),K(n)) != siteIndex(psi,n) )
          Error("MPS and MPO have different site indices in applyMPO method 'DensityMatrix'");
      }

    auto plev = 14741;

    auto res = psi;

    //Set up conjugate psi and K
    auto psic = psi;
    auto Kc = K;
    psic.dag().prime(plev);
    Kc.dag().prime(plev);

    // Make sure the original and conjugates match
    for(auto j : range1(N-1)) 
        Kc.ref(j).prime(-plev,siteIndex(Kc,psic,j));

    //Build environment tensors from the left
    if(verbose) print("Building environment tensors...");
    auto E = std::vector<ITensor>(N+1);
    E[1] = psi(1)*K(1)*Kc(1)*psic(1);
    for(int j = 2; j < N; ++j)
        {
        E[j] = E[j-1]*psi(j)*K(j)*Kc(j)*psic(j);
        //assert(order(E[j])==4);
        }
    if(verbose) println("done");

    //O is the representation of the product of K*psi in the new MPS basis
    auto O = psi(N)*K(N);

    auto rho = E[N-1] * O * dag(prime(O,plev));

    ITensor U,D;
    auto ts = tags(linkIndex(psi,N-1));
    auto spec = diagHermitian(rho,U,D,{dargs,"Tags=",ts});

    if(verbose) printfln("  j=%02d truncerr=%.2E m=%d",N-1,spec.truncerr(),dim(commonIndex(U,D)));

    res.ref(N) = dag(U);

    O = O*U*psi(N-1)*K(N-1);

    for(int j = N-1; j > 1; --j)
        {
        if(not maxdim_set)
            {
            //Infer maxdim from bond dim of original MPS
            //times bond dim of MPO
            //i.e. upper bound on order of rho
            auto cip = commonIndex(psi(j),E[j-1]);
            auto ciw = commonIndex(K(j),E[j-1]);
            auto maxdim = (cip) ? dim(cip) : 1l;
            maxdim *= (ciw) ? dim(ciw) : 1l;
            dargs.add("MaxDim",maxdim);
            }
        rho = E[j-1] * O * dag(prime(O,plev));
        ts = tags(linkIndex(psi,j-1));
        auto spec = diagHermitian(rho,U,D,{dargs,"Tags=",ts});
        O = O*U*psi(j-1)*K(j-1);
        res.ref(j) = dag(U);
        if(verbose) printfln("  j=%02d truncerr=%.2E m=%d",j,spec.truncerr(),dim(commonIndex(U,D)));
        }

    if(normalize) O /= norm(O);
    res.ref(1) = O;
    res.leftLim(0);
    res.rightLim(2);

    return res;
    }

void
fitApplyMPOImpl(Real fac,
                MPS const& x,
                MPO const& K,
                MPS& Kx,
                Sweeps const& sweeps,
                Args args)
    {
    auto N = length(x);
    auto verbose = args.getBool("Verbose",false);
    auto normalize = args.getBool("Normalize",false);

    for( auto n : range1(N) )
        if( siteIndex(Kx,n)!=siteIndex(K,x,n) )
            Error("In applyMPO with Method=Fit, guess MPS must have the same sites that the result of MPO*MPS would have");

    auto rand_plev = 43154353;
    Kx.dag().primeLinks(rand_plev);
    Kx.position(1);

    auto E = vector<ITensor>(N+2);
    E[N] = x(N)*K(N)*Kx(N);
    for(auto n = N-1; n > 2; --n)
        E[n] = E[n+1]*x(n)*K(n)*Kx(n);

    for(auto sw : range1(sweeps.nsweep()))
        {
        args.add("Sweep",sw);
        args.add("Cutoff",sweeps.cutoff(sw));
        args.add("MinDim",sweeps.mindim(sw));
        args.add("MaxDim",sweeps.maxdim(sw));
        args.add("Noise",sweeps.noise(sw));

        for(int b = 1, ha = 1; ha <= 2; sweepnext(b,ha,N))
            {
            if(verbose)
                printfln("Sweep=%d, HS=%d, Bond=(%d,%d)",sw,ha,b,b+1);

            auto lwfK = (E[b-1] ? E[b-1]*x(b) : x(b));
            lwfK *= K(b);
            auto rwfK = (E[b+2] ? E[b+2]*x(b+1) : x(b+1));
            rwfK *= K(b+1);

            auto wfK = lwfK*rwfK;
            wfK *= fac;

            if(normalize) wfK /= norm(wfK);
            auto PH = LocalOp(K(b),K(b+1),E[b-1],E[b+2]);

            wfK.dag();
            auto spec = Kx.svdBond(b,wfK,(ha==1?Fromleft:Fromright),PH,args);
 
            if(verbose)
                {
                printfln("    Trunc. err=%.1E, States kept=%s",
                         spec.truncerr(),
                         showDim(linkIndex(Kx,b)) );
                }

            if(ha == 1)
                E[b] = lwfK * Kx(b);
            else
                E[b+1] = rwfK * Kx(b+1);
            }
        }
    Kx.dag().primeLinks(-rand_plev);
    }

void
fitApplyMPOImpl(Real fac,
                MPS const& psi,
                MPO const& K,
                MPS& res,
                Args args)
    {
    if( args.defined("Maxm") )
      {
      if( args.defined("MaxDim") )
        {
        Global::warnDeprecated("Args Maxm and MaxDim are both defined. Maxm is deprecated in favor of MaxDim, MaxDim will be used.");
        }
      else
        {
        Global::warnDeprecated("Arg Maxm is deprecated in favor of MaxDim.");
        args.add("MaxDim",args.getInt("Maxm"));
        }
      }

    auto nsweep = args.getInt("Nsweep",1);
    Sweeps sweeps(nsweep);
    auto cutoff = args.getReal("Cutoff",-1);
    if(cutoff >= 0) sweeps.cutoff() = cutoff;
    auto maxdim = args.getInt("MaxDim",-1);
    if(maxdim >= 1) sweeps.maxdim() = maxdim;
    fitApplyMPOImpl(fac,psi,K,res,sweeps,args);
    }

void
fitApplyMPOImpl(MPS const& psi,
                MPO const& K,
                MPS& res,
                Args const& args)
    {
    fitApplyMPOImpl(1.,psi,K,res,args);
    }

void
applyExpH(MPS const& psi, 
          MPO const& H, 
          Real tau, 
          MPS& res, 
          Args const& args)
    {
    if(&psi == &res) Error("Must pass distinct MPS arguments to applyExpH");

    const int order = args.getInt("Order",10);

    const int N = length(res);
    const int nsweep = args.getInt("Nsweep",1);

    res.position(1);

    vector<ITensor> lastB(N+2),
                   B(N+2),
                   BH(N+2);

    B.at(N) = psi(N)*dag(prime(psi(N),"Link"));
    BH.at(N) = psi(N)*H(N)*dag(prime(psi(N)));
    for(int n = N-1; n > 2; --n)
        {
        B.at(n) = B.at(n+1)*psi(n)*dag(prime(psi(n),"Link"));
        BH.at(n) = BH.at(n+1)*psi(n)*H(n)*dag(prime(psi(n)));
        }

    lastB = B;

    MPS last(psi);

    bool up = true;

    for(int ord = order, n = 0; ord >= 1; --ord, ++n)
        {
        const Real mpofac = -tau/(1.*ord);

        if(n > 0) lastB.swap(B);

        for(int sw = 1; sw <= nsweep; ++sw)
            {
            for(int b = 1, ha = 1; ha <= 2; sweepnext(b,ha,N))
                {
                ITensor lwf,rwf,
                       lwfH,rwfH;

                if(up)
                    {
                    lwf = (B.at(b-1) ? B.at(b-1)*psi(b) : psi(b) );
                    rwf = (B.at(b+2) ? B.at(b+2)*psi(b+1) : psi(b+1));

                    lwfH = (BH.at(b-1) ? BH.at(b-1)*last(b) : last(b));
                    lwfH *= H(b);
                    rwfH = (BH.at(b+2) ? BH.at(b+2)*last(b+1) : last(b+1));
                    rwfH *= H(b+1);
                    }
                else //dn
                    {
                    lwf = (B.at(b-1) ? B.at(b-1)*dag(prime(psi(b),"Link")) : dag(prime(psi(b),"Link")));
                    rwf = (B.at(b+2) ? B.at(b+2)*dag(prime(psi(b+1),"Link")) : dag(prime(psi(b+1),"Link")));

                    lwfH = (BH.at(b-1) ? BH.at(b-1)*dag(prime(last(b))) : dag(prime(last(b))));
                    lwfH *= H(b);
                    rwfH = (BH.at(b+2) ? BH.at(b+2)*dag(prime(last(b+1))) : dag(prime(last(b+1))));
                    rwfH *= H(b+1);
                    }

                auto wf = noPrime(lwf*rwf) + mpofac*noPrime(lwfH*rwfH);
                if(!up) wf.dag();

                res.svdBond(b,wf,(ha==1?Fromleft:Fromright),args+Args("UseSVD",true));

                if(up)
                    {
                    if(ha == 1)
                        {
                        B.at(b) = lwf * dag(prime(res(b),"Link"));
                        BH.at(b) = lwfH * dag(prime(res(b)));
                        }
                    else
                        {
                        B.at(b+1) = rwf * dag(prime(res(b+1),"Link"));
                        BH.at(b+1) = rwfH * dag(prime(res(b+1)));
                        }
                    }
                else //dn
                    {
                    if(ha == 1)
                        {
                        B.at(b) = lwf * res(b);
                        BH.at(b) = lwfH * res(b);
                        }
                    else
                        {
                        B.at(b+1) = rwf * res(b+1);
                        BH.at(b+1) = rwfH * res(b+1);
                        }
                    }
                }
            }

        last = res;

        up = !up;

        } // for ord

    }

//
// Deprecated
//

//
// For now this is unsupported
// We can consider bringing it back for example
// as a way to get a default starting MPS
// for fitApplyMPOImpl
//

//void 
//zipUpApplyMPOImpl(MPS const& psi, 
//                  MPO const& K, 
//                  MPS& res, 
//                  Args const& args)
//    {
//    const
//    bool allow_arb_position = args.getBool("AllowArbPosition",false);
//
//    if(&psi == &res)
//        Error("psi and res must be different MPS instances");
//
//    auto N = length(psi);
//    if(length(K) != N) 
//        Error("Mismatched N in ApplyMPO() ZipUp method");
//
//    if(!itensor::isOrtho(psi) || itensor::orthoCenter(psi) != 1)
//        Error("Ortho center of psi must be site 1");
//
//    if(!allow_arb_position && (!itensor::isOrtho(K) || itensor::orthoCenter(K) != 1))
//        Error("Ortho center of K must be site 1");
//
//#ifdef DEBUG
//    checkQNs(psi);
//    checkQNs(K);
//    /*
//    cout << "Checking divergence in zip" << endl;
//    for(int i = 1; i <= N; i++)
//	div(psi(i));
//    for(int i = 1; i <= N; i++)
//	div(K(i));
//    cout << "Done Checking divergence in zip" << endl;
//    */
//#endif
//
//    res = psi; 
//    res.replaceTags("0","4","Link");
//    res.replaceTags("0","1","Site");
//
//    ITensor clust,nfork;
//    vector<int> midsize(N);
//    int maxdim = 1;
//    for(int i = 1; i < N; i++)
//        {
//        if(i == 1) { clust = psi(i) * K(i); }
//        else { clust = nfork * (psi(i) * K(i)); }
//        if(i == N-1) break; //No need to SVD for i == N-1
//
//        Index oldmid = linkIndex(res,i); assert(oldmid.dir() == Out);
//        nfork = ITensor(linkIndex(psi,i),linkIndex(K,i-1),oldmid);
//        //if(clust.iten_size() == 0)	// this product gives 0 !!
//	    //throw ResultIsZero("clust.iten size == 0");
//        denmatDecomp(clust, res.ref(i), nfork,Fromleft,args);
//        Index mid = commonIndex(res(i),nfork);
//        //assert(mid.dir() == In);
//        mid.dag();
//        midsize[i] = dim(mid);
//        maxdim = std::max(midsize[i],maxdim);
//        assert(linkIndex(res,i+1).dir() == Out);
//        res.ref(i+1) = ITensor(mid,prime(res.sites()(i+1)),linkIndex(res,i+1));
//        }
//    nfork = clust * psi(N) * K(N);
//    //if(nfork.iten_size() == 0)	// this product gives 0 !!
//    //throw ResultIsZero("nfork.iten size == 0");
//
//    res.svdBond(N-1,nfork,Fromright,args);
//    res.noPrime("Link");
//    res.replaceTags("1","0","Site");
//    res.position(1);
//    } //void zipUpApplyMPOImpl

//
// These versions calculate |res> = |psiA> + mpofac*H*|psiB>
// Currently they are unsupported

//Real
//fitApplyMPOImpl(MPS const& psiA, 
//                Real mpofac,
//                MPS const& psiB,
//                MPO const& K,
//                MPS& res,
//                Args const& args)
//    {
//    return fitApplyMPOImpl(1.,psiA,mpofac,psiB,K,res,args);
//    }
//
//
//Real
//fitApplyMPOImpl(Real mpsfac,
//                MPS const& psiA, 
//                Real mpofac,
//                MPS const& psiB,
//                MPO const& K,
//                MPS& res,
//                Args const& args)
//    {
//    if(&psiA == &res || &psiB == &res)
//        {
//        Error("fitApplyMPOImpl: Result MPS cannot be same as an input MPS");
//        }
//    auto N = length(psiA);
//    auto nsweep = args.getInt("Nsweep",1);
//
//    res.position(1);
//
//    vector<ITensor> B(N+2),
//                   E(N+2);
//
//    B.at(N) = psiA(N)*dag(prime(res(N),"Link"));
//    E.at(N) = psiB(N)*K(N)*dag(prime(res(N)));
//    for(int n = N-1; n > 2; --n)
//        {
//        B.at(n) = B.at(n+1)*psiA(n)*dag(prime(res(n),"Link"));
//        E.at(n) = E.at(n+1)*psiB(n)*K(n)*dag(prime(res(n)));
//        }
//
//
//    for(int sw = 1; sw <= nsweep; ++sw)
//        {
//        for(int b = 1, ha = 1; ha <= 2; sweepnext(b,ha,N))
//            {
//            ITensor lwf = (B.at(b-1) ? B.at(b-1)*psiA(b) : psiA(b));
//            ITensor rwf = (B.at(b+2) ? psiA(b+1)*B.at(b+2) : psiA(b+1));
//
//            ITensor lwfK = (E.at(b-1) ? E.at(b-1)*psiB(b) : psiB(b));
//            lwfK *= K(b);
//            ITensor rwfK = (E.at(b+2) ? E.at(b+2)*psiB(b+1) : psiB(b+1));
//            rwfK *= K(b+1);
//
//            ITensor wf = mpsfac*noPrime(lwf*rwf) + mpofac*noPrime(lwfK*rwfK);
//            wf.noPrime();
//
//            res.svdBond(b,wf,(ha==1?Fromleft:Fromright),args+Args("UseSVD",true));
//
//            if(ha == 1)
//                {
//                B.at(b) = lwf * dag(prime(res(b),"Link"));
//                E.at(b) = lwfK * dag(prime(res(b)));
//                }
//            else
//                {
//                B.at(b+1) = rwf * dag(prime(res(b+1),"Link"));
//                E.at(b+1) = rwfK * dag(prime(res(b+1)));
//                }
//            }
//        }
//
//    auto olp = B.at(3);
//    olp *= psiA(2);
//    olp *= dag(prime(res(2),"Link"));
//    olp *= psiA(1);
//    olp *= dag(prime(res(1),"Link"));
//
//    return olp.elt();
//    }

} //namespace itensor
