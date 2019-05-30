#include "itensor/all.h"

using namespace itensor;

ITensor
makeSp(Index const& s)
    {
    auto Sp = ITensor(s,prime(s));
    Sp.set(s=2,prime(s)=1, 1);
    return Sp;
    }

ITensor
makeSm(Index const& s)
    {
    auto Sm = ITensor(s,prime(s));
    Sm.set(s=1,prime(s)=2,1);
    return Sm;
    }

ITensor
makeSz(Index const& s)
    {
    auto Sz = ITensor(s,prime(s));
    Sz.set(s=1,prime(s)=1, 0.5);
    Sz.set(s=2,prime(s)=2,-0.5);
    return Sz;
    }

int main()
    {
    //
    // Two-site wavefunction
    // initialized to a singlet
    //
    
    auto s1 = Index(2, "s1");
    auto s2 = Index(2, "s2");

    auto psi = ITensor(s1,s2); //default initialized to zero

    psi.set(s1=1,s2=2, 1./sqrt(2));
    psi.set(s1=2,s2=1,-1./sqrt(2));

    PrintData(psi);

    //EXIT //uncomment to exit here

    //
    // Single-site operators
    //

    auto Sz1 = makeSz(s1);
    auto Sz2 = makeSz(s2);
    auto Sp1 = makeSp(s1);
    auto Sp2 = makeSp(s2);
    auto Sm1 = makeSm(s1);
    auto Sm2 = makeSm(s2);

    PrintData(Sz1);
    PrintData(Sp1);
    PrintData(Sm1);

    //EXIT //uncomment to exit here

    //
    // Two-site Heisenberg Hamiltonian
    //

    auto H = Sz1*Sz2 + 0.5*Sp1*Sm2 + 0.5*Sm1*Sp2;

    //
    // Energy expectation value
    //

    auto psidag = dag(prime(psi));
    auto E = elt(psidag * H * psi);

    Print(E);

    return 0;
    }


// psi = 
// ITensor ord=2: (2|id=862|s1) (2|id=140|s2) 
// {norm=1.00 (Dense Real)}
// (2,1) -0.707107
// (1,2) 0.7071068

// Sz1 = 
// ITensor ord=2: (2|id=862|s1) (2|id=862|s1)' 
// {norm=0.71 (Dense Real)}
// (1,1) 0.5000000
// (2,2) -0.500000

// Sp1 = 
// ITensor ord=2: (2|id=862|s1) (2|id=862|s1)' 
// {norm=1.00 (Dense Real)}
// (2,1) 1.0000000

// Sm1 = 
// ITensor ord=2: (2|id=862|s1) (2|id=862|s1)' 
// {norm=1.00 (Dense Real)}
// (1,2) 1.0000000

// E = -0.75

