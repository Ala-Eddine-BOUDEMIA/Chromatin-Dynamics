import plotly.express as px

def histogram(file, title, *args, **kwargs):

	fig = px.histogram(
		file, x=kwargs.get('x',None), 
		y=kwargs.get('y',None), title=title)
	
	fig.show()

def boxplot(file, title, *args, **kwargs):

	fig = px.box(
		file, x=kwargs.get('x',None), 
		y=kwargs.get('y',None), title=title)
	
	fig.show()