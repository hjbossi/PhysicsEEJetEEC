import ROOT, yaml, math

with open('settings.yaml', 'r') as stream:
    try:
        setup = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)
        exit()

NominalFile = ROOT.TFile(setup['nominal']['file'])
NominalHistogram = NominalFile.Get(setup['nominal']['histogram'])
NominalHistogram.SetName("NominalHistogram")

OutputFile = ROOT.TFile(setup['output']['file'], 'recreate')
TotalHistogram = NominalHistogram.Clone()
TotalHistogram.SetNameTitle(setup['output']['histogram'] + '_Total', f"Total systematics for {setup['output']['histogram']}")
TotalHistogram.Reset()

# Loop over systematics sources
for Tag in setup['systematics']:
    print(f".....Processing the systematic of type = {Tag}")
    Item = setup['systematics'][Tag]

    # Get the base histogram
    if Item['base_type'].casefold() == 'Nominal'.casefold():
        Base = NominalHistogram.Clone(setup['output']['histogram'] + f"_{Tag}_Base")
    elif Item['base_type'].casefold() == 'File'.casefold():
        BaseFile = ROOT.TFile(Item['base_file'])
        OutputFile.cd()
        Base = BaseFile.Get(Item['base_histogram']).Clone(setup['output']['histogram'] + f"_{Tag}_Base")
        BaseFile.Close()
    else:
        print(f"Base type {Item['base_type']} not found!")
    Base.SetNameTitle(setup['output']['histogram'] + f"_{Tag}_Base", f"Base for {Tag}")

    # Get the variant histogram
    if Item['variant_type'].casefold() == 'Nominal'.casefold():
        OutputFile.cd()
        Variant = NominalHistogram.Clone(setup['output']['histogram'] + f"_{Tag}_Variant")
    elif Item['variant_type'].casefold() == 'File'.casefold():
        VariantFile = ROOT.TFile(Item['variant_file'])
        OutputFile.cd()
        Variant = VariantFile.Get(Item['variant_histogram']).Clone(setup['output']['histogram'] + f"_{Tag}_Variant")
        VariantFile.Close()
    elif Item['variant_type'].casefold() == 'Scale'.casefold():
        print("HERE")
        OutputFile.cd()
        Variant = NominalHistogram.Clone(setup['output']['histogram'] + f"_{Tag}_Variant")
        Variant.Scale(Item['variant_scale'])
    elif Item['variant_type'].casefold() == '1DRatio'.casefold():
        OutputFile.cd()
        Variant = NominalHistogram.Clone(setup['output']['histogram'] + f"_{Tag}_Variant")
        VariantFile = ROOT.TFile(Item['variant_file'])
        Numerator = VariantFile.Get(Item['variant_numerator'])
        Denominator = VariantFile.Get(Item['variant_denominator'])
        Scale = Item['variant_scale'] if 'variant_scale' in Item else 1.00
        for iX in range(0, NBinsX + 2):   # include also underflow and overflow (0 and N+1)
            for iY in range(0, NBinsY + 2):
                N = Numerator.GetBinContent(iX)
                D = Denominator.GetBinContent(iX)
                if D != 0:
                    Correction = N / D
                else:
                    Correction = 1
                if Correction <= 0:
                    Correction = 1
                Variant.SetBinContent(iX, iY, Variant.GetBinContent(iX, iY) * ((Correction - 1) * Scale + 1))
        VariantFile.Close()
    else:
        print(f"Base type {Item['base_type']} not found!")
    Variant.SetNameTitle(setup['output']['histogram'] + f"_{Tag}_Variant", f"Variant for {Tag}")

    # Now we add whatever we want into the total histogram
    NBinsX = TotalHistogram.GetNbinsX()
    NBinsY = TotalHistogram.GetNbinsY()
    for iX in range(0, NBinsX + 2):   # include also underflow and overflow (0 and N+1)
        for iY in range(0, NBinsY + 2):
            Current = TotalHistogram.GetBinContent(iX, iY)
            Extra = Variant.GetBinContent(iX, iY) - Base.GetBinContent(iX, iY)
            if(Base.GetBinContent(iX, iY) != 0):
                RelError = Variant.GetBinContent(iX, iY)/ Base.GetBinContent(iX, iY)
                print(f"Variant bin content is {Variant.GetBinContent(iX, iY)}, Base bin content is {Base.GetBinContent(iX, iY)}, Percent error of bin {iX},{iY} variation is {RelError}")
            Current = math.sqrt(Current * Current + Extra * Extra)
            TotalHistogram.SetBinContent(iX, iY, Current)

    # Write the individual ones to file
    OutputFile.cd()
    Base.Write()
    Variant.Write()

# Write things out
OutputFile.cd()
TotalHistogram.Write()
NominalHistogram.Write()

print('Done!')

OutputFile.Close()

print('Done!')

# NominalFile.Close()




