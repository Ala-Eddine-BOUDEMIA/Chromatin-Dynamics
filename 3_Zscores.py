import pandas as pd
import plotly.express as px

from sklearn.preprocessing import StandardScaler

import Tools
import Config

def z_scores(
	rand, top88, cv_counts, raw_counts, top1000,
	counts_norm, tissue_counts, normal, filtered_counts,
	nrcv_counts, Icv, Iraw, Irand, Inrcv, Itop88, Itop1000,
	Inormal, Itissue, Inormalized, Ifiltered, IwoTissues, Pcv,
	Praw, Prand, Pnrcv, Ptop88, Ptop1000, Pnormal, Ptissue, Pnormalized, 
	Pfiltered, PwoTissues, counts_without_tissues):
		
	counts = [raw_counts, filtered_counts, counts_norm, normal, 
        counts_without_tissues, top1000, top88, cv_counts, nrcv_counts]
		
	tissue_files = sorted([f for f in tissue_counts.iterdir() if f.is_file()])
	for path in tissue_files:
		counts.append(path)
		
	rand_files = sorted([f for f in rand.iterdir() if f.is_file()])
	for path in rand_files:
		counts.append(path)
	
	images = [Iraw, Ifiltered, Inormalized, Inormal, IwoTissues, 
		Itop1000, Itop88, Icv, Inrcv]
    
	htmls = [Praw, Pfiltered, Pnormalized, Pnormal, PwoTissues,
		Ptop1000, Ptop88, Pcv, Pnrcv]

	for i in range(len(tissue_files)):
		tissue_name = str(tissue_files[i]).split("/")[-1].split(".")[0]
		link_img = Itissue.joinpath(tissue_name)
		link_p = Ptissue.joinpath(tissue_name)
		
		Tools.create_folder(link_img)
		Tools.create_folder(link_p)
		
		images.append(link_img.joinpath(tissue_name + ".png"))
		htmls.append(link_p.joinpath(tissue_name + ".html"))
		
	for i in range(len(rand_files)):
		link_img = Irand.joinpath("random" + str(i))
		link_html = Prand.joinpath("random" + str(i))
		
		Tools.create_folder(link_img)
		Tools.create_folder(link_html)
		
		images.append(link_img.joinpath("random" + str(i) + ".png"))
		htmls.append(link_html.joinpath("random" + str(i) + ".html"))

	for count, i, h in zip(counts, images, htmls):
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

if __name__ == '__main__':	
	z_scores(
		rand = Config.args.rand,
		top88 = Config.args.top88, 
		cv_counts = Config.args.cv,
		raw_counts = Config.args.bf, 
		top1000 = Config.args.top1000,
		counts_norm = Config.args.norm, 
		tissue_counts = Config.args.tissue,
		normal = Config.args.onlyNormal,
		filtered_counts = Config.args.af,
		nrcv_counts = Config.args.nonRcv, 
		Icv = Config.args.IzscoreCV,
		Iraw = Config.args.IzscoreRaw,
		Irand = Config.args.IzscoreRand,
		Inrcv = Config.args.IzscoreNrCV,
		Itop88 = Config.args.IzscoreTop88,
		Itop1000 = Config.args.IzscoreTop,
		Inormal = Config.args.IzscoreNormal,
		Itissue = Config.args.IzscoreTissue,
		Inormalized = Config.args.IzscoreNorm,
		Ifiltered = Config.args.IzscoreFiltered,
		IwoTissues = Config.args.IzscoreWoTissues,
		Pcv = Config.args.PzscoreCV,
		Praw = Config.args.PzscoreRaw,
		Prand = Config.args.PzscoreRand,
		Pnrcv = Config.args.PzscoreNrCV,
		Ptop88 = Config.args.PzscoreTop88,
		Ptop1000 = Config.args.PzscoreTop,
		Pnormal = Config.args.PzscoreNormal,
		Ptissue = Config.args.PzscoreTissue,
		Pnormalized = Config.args.PzscoreNorm,
		Pfiltered = Config.args.PzscoreFiltered,
		PwoTissues = Config.args.PzscoreWoTissues,
		counts_without_tissues = Config.args.WoTissues)