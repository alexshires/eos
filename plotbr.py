#! /usr/bin/env python



import array

def get_br_from_file( filename, xv, xe, yc, yd, yu ) :
    f = open( filename )
    i = 0
    for line in f :
        i += 1
        if i < 3:
            continue
        ls = line.split()
        xmin = float(ls[0])
        xmax = float(ls[1])
        xv.append(xmin+0.5*(xmax-xmin))
        xe.append(0.5*(xmax-xmin))
        yc.append(float(ls[2]))
        yd.append(float(ls[5]))
        yu.append(float(ls[6]))
    f.close()


def get_dbr_from_file( filename, xv, xe, yc, yd, yu ) :
    f = open( filename )
    i = 0
    for line in f :
        i += 1
        if i < 3:
            continue
        ls = line.split()
        xv.append(float(ls[0]))
        xe.append(float(0))
        yc.append(float(ls[1]))
        yd.append(float(ls[4]))
        yu.append(float(ls[5]))
    f.close()


qsqvals = []
qsqerrs = []
pimmbrcvals = []
pimmbruerrs = []
pimmbrderrs = []




get_dbr_from_file( 'pimumu.large.data' 
        , qsqvals, qsqerrs
        , pimmbrcvals, pimmbruerrs, pimmbrderrs 
        )
get_dbr_from_file( 'pimumu.low.data' 
        , qsqvals, qsqerrs
        , pimmbrcvals, pimmbruerrs, pimmbrderrs 
        )


kmmbrcvals = []
kmmbruerrs = []
kmmbrderrs = []

get_dbr_from_file( 'kmumu.large.data' 
        , [], []
        , kmmbrcvals, kmmbruerrs, kmmbrderrs 
        )
get_dbr_from_file( 'kmumu.low.data' 
        , [], []
        , kmmbrcvals, kmmbruerrs, kmmbrderrs 
        )


qsqbins = []
qsqends = []
pimmbins = []
pimmbinup = []
pimmbindown = []
kmmbins = []
kmmbinup = []
kmmbindown = []

get_br_from_file( 'pimumu.large.value.data' 
        , qsqbins, qsqends
        , pimmbins, pimmbinup, pimmbindown 
        )
get_br_from_file( 'pimumu.low.value.data' 
        , qsqbins, qsqends
        , pimmbins, pimmbinup, pimmbindown
        )
get_br_from_file( 'kmumu.large.value.data' 
        , [], []
        , kmmbins, kmmbinup, kmmbindown 
        )
get_br_from_file( 'kmumu.low.value.data' 
        , [], []
        , kmmbins, kmmbinup, kmmbindown
        )






import numpy as np

qsqvals = np.array( qsqvals )
qsqerrs = np.array( qsqerrs )
pimmbrcvals = np.array (pimmbrcvals)
pimmbrderrs = np.array (pimmbrderrs)
pimmbruerrs = np.array (pimmbruerrs)
kmmbrcvals = np.array (kmmbrcvals)
kmmbrderrs = np.array (kmmbrderrs)
kmmbruerrs = np.array (kmmbruerrs)


qsqbins = np.array( qsqbins ) 
qsqends = np.array( qsqends ) 
pimmbins = np.array( pimmbins ) 
pimmbinup = np.array( pimmbinup ) 
pimmbindown = np.array( pimmbindown ) 
kmmbins = np.array( kmmbins ) 
kmmbinup = np.array( kmmbinup ) 
kmmbindown = np.array( kmmbindown ) 

print qsqbins, pimmbins, kmmbins


from matplotlib import rc, rc_file
rc_file('/home/alexshires/.config/matplotlib/matplotlibrc')
from matplotlib import rc
rc('text', usetex=True)

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages("hmmbr.pdf")

plt.semilogy()

kmm = plt.errorbar( qsqbins, kmmbins
                   , xerr=qsqends, yerr=[kmmbindown,kmmbinup], fmt=' '
                   , label="$K^{+}\mu^{+}\mu^{-}$", color='green'  
                   )

pmm = plt.errorbar( qsqbins, pimmbins
             , xerr=qsqends, yerr=[pimmbindown, pimmbinup], fmt=' '
             , label="$\pi^{+}\mu^{+}\mu^{-}$", color='blue'
              )

print kmm, pmm

plt.ylabel('BR($B^{+} \\rightarrow h^{+}\mu^{+}\mu^{-}$)', ha='right', y=1)
plt.xlabel('$q^{2}$ [GeV$^{2}$]', ha='right', x=1)
plt.legend()
plt.minorticks_on()
pp.savefig()
#calculate ratio


ratiobins = pimmbins / kmmbins / (0.216**2)
a = kmmbinup / kmmbins
aa = a * a
b = pimmbinup / pimmbins
bb = b * b
ratiobinup = ratiobins * np.sqrt(aa + bb) 
a = kmmbindown / kmmbins
aa = a * a
b = pimmbindown / pimmbins
bb = b * b
ratiobindown = ratiobins * np.sqrt(aa + bb)

print ratiobins[0], ratiobinup[0], ratiobindown[0]
print ratiobins[1], ratiobinup[1], ratiobindown[1]




plt.clf()
plt.semilogy()
plt.fill_between(qsqvals, pimmbrcvals-pimmbrderrs, pimmbrcvals+pimmbruerrs,
                 facecolor='blue', interpolate=True)
plt.plot(qsqvals, pimmbrcvals, label="$\\pi^{+}\mu^{+}\mu^{-}$",
         color='black')

plt.ylabel('BR($B^{+} \\rightarrow h^{+}\\mu^{+}\\mu^{-}$)', ha='right', y=1)
plt.xlabel('$q^{2}$ [GeV$^{2}$]', ha='right', x=1)
plt.legend()
plt.minorticks_on()
pp.savefig()
plt.fill_between(qsqvals, kmmbrcvals-kmmbrderrs, kmmbrcvals+kmmbruerrs,
                 facecolor='green', interpolate=True)
plt.plot(qsqvals, kmmbrcvals, label="$K^{+}\mu^{+}\mu^{-}$",
         color='black', linestyle='-.')



plt.ylabel('BR($B^{+} \\rightarrow h^{+}\\mu^{+}\\mu^{-}$)', ha='right', y=1)
plt.xlabel('$q^{2}$ [GeV$^{2}$]', ha='right', x=1)
plt.legend()
plt.minorticks_on()
pp.savefig()



plt.clf()

ratiocvals = kmmbrcvals / pimmbrcvals

a = kmmbruerrs / kmmbrcvals
aa = a * a
b = pimmbruerrs / pimmbrcvals
bb = b * b

ratiouperr = ratiocvals + ratiocvals * np.sqrt(aa + bb) 
ratiodoerr = ratiocvals - ratiocvals * np.sqrt(aa + bb) 

plt.fill_between(qsqvals, ratiouperr, ratiodoerr, 
                 facecolor='gray', interpolate=True)
plt.plot(qsqvals, ratiocvals, label = "ratio", color='black', linestyle='-.')
plt.ylabel("$B(K^{+}\mu^{+}\mu^{-}) / B(\pi^{+}\mu^{+}\mu^{-})$", ha='right', y=1)
plt.xlabel('$q^{2}$ [GeV$^{2}$]', ha='right', x=1)
#plt.ylim([0,3])
plt.legend()
plt.minorticks_on()
pp.savefig()

plt.clf()

nockmuperr = ratiouperr * (0.216 ** 2)
nockmdoerr = ratiodoerr * (0.216 ** 2)
nockmcvals = ratiocvals * (0.216 ** 2)
plt.fill_between(qsqvals, nockmuperr, nockmdoerr, 
                 facecolor='gray', interpolate=True)
plt.plot(qsqvals, nockmcvals, label = "ratio * $|V_{td}/V_{ts}|^{2}$", color='blue')
plt.ylabel("$F_K(q^{2}) / F_{\pi}(q^{2})$", ha='right', y=1)
plt.xlabel('$q^{2}$ [GeV$^{2}$]', ha='right', x=1)
#plt.ylim([0,0.5])
plt.legend()
plt.minorticks_on()
pp.savefig()

pp.close()





