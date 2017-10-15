using DFControl

bands = read_qe_bands_file("/home/ponet/Documents/PhD/GeTe/SOC/GeTe_bands.out")
plot(plot_qe_bands("/home/ponet/Documents/PhD/GeTe/SOC/GeTe_bands.out"),plot_qe_kpdos("/home/ponet/Documents/PhD/GeTe/NSOC/pwo.pdos_atm#2(Ge)_wfc#2(p)"))
