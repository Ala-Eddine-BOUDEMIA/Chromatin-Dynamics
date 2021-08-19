import gc
import pandas as pd
import plotly.express as px

from sklearn.preprocessing import StandardScaler

import Tools
import Config

def z_scores(
	counts, rand, tissues,
	Icounts, Irand, Itissue,
	Pcounts, Prand, Ptissue):

	tissue_files = sorted([f for f in tissues.glob("**/*.tsv") if f.is_file()])
	for path in tissue_files:
		counts.append(path)
	
	"""rand_files = sorted([f for f in rand.glob("**/*.tsv") if f.is_file()])
	for path in rand_files:
		counts.append(path)"""

	for i in range(len(tissue_files)):
		tissue_name = str(tissue_files[i]).split("/")[-1].split(".")[0]
		link_img = Itissue.joinpath(tissue_name)
		link_p = Ptissue.joinpath(tissue_name)
		
		Tools.create_folder(link_img)
		Tools.create_folder(link_p)
		
		Icounts.append(link_img.joinpath(tissue_name + ".png"))
		Pcounts.append(link_p.joinpath(tissue_name + ".html"))

	"""for i in range(len(rand_files)):
		link_img = Irand.joinpath("random" + str(i))
		link_html = Prand.joinpath("random" + str(i))
		
		Tools.create_folder(link_img)
		Tools.create_folder(link_html)
		
		Icounts.append(link_img.joinpath("random" + str(i) + ".png"))
		Pcounts.append(link_html.joinpath("random" + str(i) + ".html"))"""

	for count, i, h in zip(counts, Icounts, Pcounts):
		c = pd.read_csv(count, header = 0, index_col = 0, sep = "\t")
		 
		scaler = StandardScaler()
		std_counts = pd.DataFrame(scaler.fit_transform(c), 
			index = c.index, columns = c.columns)

		df = pd.DataFrame(index = c.index)
		mean = pd.DataFrame(std_counts.mean(axis = 1), columns = ["mean"])
		std = pd.DataFrame(std_counts.std(axis = 1), columns = ['std']) 

		df = df.join(mean)
		df = df.join(std)
		df = df.sort_values("mean")

		fig = px.scatter(
			df, x = "std", y = "mean", 
			hover_data = [df.index], 
			title = "z_scores")
		
		fig.write_html(str(h))
		fig.write_image(str(i), width = 2048, height = 1024)
		#fig.show()

		del(c)
		del(df)
		del(fig)
		gc.collect()

if __name__ == '__main__':	
	z_scores(
		counts = Config.counts,
		rand = Config.args.rand,
		tissues = Config.args.tissue,
		Icounts = Config.zscores_imgs,
		Irand = Config.args.IzscoreRand,
		Itissue = Config.args.IzscoreTissue,
		Pcounts = Config.zscores_htmls,
		Prand = Config.args.PzscoreRand,
		Ptissue = Config.args.PzscoreTissue)