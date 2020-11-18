import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt

def plotAD(dat, out, stat='frequency', binwidth=0.1, y=None, x=None):
	sns.set_theme(style="ticks")
	f, ax = plt.subplots(figsize=(7, 5))
	sns.despine()
	
	sns_plot = sns.histplot(
		dat,
		x="AD", 
		multiple="stack",
		palette="colorblind",
		edgecolor=".3",
		linewidth=.5,
		stat=stat,
		binwidth=binwidth
	)
	plt.ylim(0, y)
	plt.xlim(0, x)
	oname=out + "_hist.pdf"
	#plt.show()
	sns_plot.figure.savefig(oname, dpi=400)

def plotADclassified(dat, out, cls, stat='frequency', binwidth=0.1, y=None, x=None):
	sns.set_theme(style="ticks")
	f, ax = plt.subplots(figsize=(7, 5))
	sns.despine()
	
	sns_plot = sns.histplot(
		dat,
		x="AD", hue=cls,
		multiple="stack",
		palette="colorblind",
		edgecolor=".3",
		linewidth=.5,
		stat=stat,
		binwidth=binwidth
	)
	plt.ylim(0, y)
	plt.xlim(0, x)

	
	oname=out + "_hist.pdf"
	#plt.show()
	sns_plot.figure.savefig(oname, dpi=400)