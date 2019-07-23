<!DOCTYPE html>
<html lang="en">

	<head>

		<meta http-equiv="content-type" content="text/html; charset=UTF-8">
		<title>reStrainingOrder Summary Report - {{filename}}</title>

		<style>
                  body {
                          font-family: Arial, sans-serif;
                          font-size:14px;
                          padding:0 20px 20px;
						  height: 100%;
                  }
				  html{
					  height: 100%;
					  margin: 0px;
				  }
				  .container {
                          margin:0 auto;
                          max-width:1200px;
                  }                  
                  .header img {
                          float:left;
                  }
                  .header h1 {
                          margin: 20px 0 10px;
                  }
                  .header img {
                          padding: 0 20px 20px 0;
                  }
                  .subtitle {
                          margin-top:120px;
                          float:right;
                          text-align:right;
                  }
                  .header_subtitle h3,
                  .header_subtitle p {
                          margin:0;
                  }
                  h1 {
                          font-size: 3.2em;
                  }
                  h2 {
                          font-size:2.2em;
                  }
                  h3 {
                          font-size:1.4em;
                  }
                  h2, h3, hr {
                          clear:both;
                  }
                  hr {
                        border-top:1px solid #CCC;
                        border-bottom:1px solid #F3F3F3;
                        border-left:0;
                        border-right:0;
                        height:0;
                  }
				  .bargraph {
						width: 100%;
						height: 100%;						
				  }
				  .plotly_table {
                         float:right;
                         width:600px;
                         max-width:100%;
                    }
                  .data {
                          float:left;
                          width:500px;
                          max-width:100%;
                          margin-right:30px;
                          border:1px solid #CCC;
                          border-collapse:separate;
                          border-spacing: 0;
                          border-left:0;
                          -webkit-border-radius:4px;
                          -moz-border-radius:4px;
                          border-radius:4px;
                  }
                  .data th, .data td {
                          border-left:1px solid #CCC;
                          border-top:1px solid #CCC;
                          padding: 5px 7px;
                  }
                  .data tr:first-child th,
                  .data tr:first-child td {
                          border-top:0;
                  }
                  .data tr:last-child th,
                  .data tr:last-child td {
                          border-bottom: 2px solid #666;
                  }
                  .plot {
                          width:650px;
                          max-width:100%;
                          float:left;
                          margin-bottom:30px;
                  }

                  .fullWidth_plot {
                          height: 600px;
                  }

                  .data th {
                          text-align:left;
                  }
                  .data td {
                          text-align:right;
                  }
                footer {
                    color:#999;
                }
                footer a {
					color:#999;
                }
                .error-msg {
                    color: #a94442;
					background-color: #f2dede;
					border: 1px solid #ebccd1;
					padding: 15px;
					margin-bottom: 20px;
					border-radius: 4px;
                }
				.error-msg h3  { margin: 0; }
				.error-msg pre { margin: 0; }
          </style>

	<!-- Plotly.js -->
	{{plotly_goes_here}}
 	This will need to be replaced by the plot.ly library itself
 	{{plotly_goes_here}}
  
	
	</head>


	<body>
	<div class="container">
		<div class="header">
		 
			<h1>reStrainingOrder Summary Report</h1>
	
			<div class="subtitle">
				<h3>{{filename}}</h3>
				<p>Data processed at {{time}} on {{date}}</p>
			</div>
			
		</div>
	
		<hr id="header_hr">
	
	

	<h2>Strain Compatibility Scores</h2>
		<table>
			<tbody>
				<tr>
					<td style="vertical-align:top">
						<table class="data">
							{{strain_compatibility_content}}
						</table>
					</td>
					<td>
						<div id="Strain_Scores_Report"><!-- Plotly chart will be drawn inside this DIV --> </div>
					</td>						
				</tr>
			</tbody>
		</table>

			
	<hr>

	<h2>Potential Hybrid Strains - Compatibility Scores</h2>
	<table>
		<tbody>
			<tr>
				<td style="vertical-align:top">
					<table class="data">
						{{hybrid_compatibility_content}}
					</table>			
				</td>
				<td>
					<div id="Hybrid_Allelic_Ratio_Plot" class="bargraph"> <!--Plotly chart will be drawn inside this DIV --> </div>
				</td>
			</tr>
		</tbody>
	</table>

	<hr>

	

	<!-- This section is optional -->

	
	<!-- ######################### PLOT.LY plotting code below ################################################################# -->

	<!-- Strain Scores Plot -->
	<script>
		var data = [
		  {
			y: [{{strain_scores_strains_plotly}}],
			x: [{{strain_scores_percentages_plotly}}],
			type: 'bar',
			marker: {
				color: '#0d233a',
				<!-- color: ['#0d233a', '#2f7ed8','#2f7ed8'], -->
	 			<!-- colors: ['#0d233a', '#2f7ed8','#8bbc21','#910000','#1aadce','#492970','#f28f43','#77a1e5','#c42525','#a6c96a'], -->
   				line: {
     				color: 'black',
      				width: 1
    			},
   			},
			orientation: 'h'
			}
		];
		
		var layout = {
			<!-- paper_bgcolor: '#7f7f7f', -->
  			<!-- plot_bgcolor: '#c7c7c7',  -->
  			height: 1000,
			width: 700,
			shapes: [{
				type: 'line',
    			x0: '100',
    			y0: 0,
				x1: '100',
				yref: 'paper',
				y1: 1,
				line: {
					color: 'grey',
					width: 1.5,
					dash: 'dot'
				}
			}],
  			margin: {
			    l: 250,
			    r: 50,
			    b: 50,
			    t: 0,
			    pad: 0,
			},
			xaxis: {
				automargin: true,
				range: [0, 100],
				title:'Strain Compatibility [in %]',
				titlefont: {
	  				size: 18, 
	  			},
			},
			yaxis:{
				autorange: "reversed",
				automargin: true,
				ticks: 'outside',
				tickcolor: 'rgba(0,0,0,0)',
			}
		
		};
		var optionsStrains= {
                        displaylogo: false,
						responsive: true,
						modeBarButtonsToRemove:
							['zoom2d',
							'pan',
							'pan2d',
							'resetScale2d',
							'hoverClosestCartesian',
							'hoverCompareCartesian',
							'toggleSpikelines']
                        ,
                        toImageButtonOptions: {
							filename: 'reStrainingOrder Strain Summary',
							width: 800,  <!-- width: stacksDivPercentage._fullLayout.width, height: stacksDivPercentage._fullLayout.height -->
							height: 600,
							format: 'png'
                        }
                };

		Plotly.newPlot('Strain_Scores_Report', data, layout, optionsStrains);
		
	</script>




	
	<!-- Hybrid Scores Plot including Allelic Ratios-->
	<script>
	  	var trace1 = {
	  		x: [{{hybrid_percentages_strain1_plotly}}], 
	  		y: [{{hybrid_compatibility_strains_plotly}}],
			  
			
	  		name: 'Strain1', 
			marker: {
				color: '#2f7ed8',
				line:{
					width: 1,
					color:'black',
				},
			}, 
		  	type: 'bar',
		  	orientation: 'h'
	  	};
	  
		 
		var trace2 = {
	  		x: [{{hybrid_percentages_strain2_plotly}}], 
	  		y: [{{hybrid_compatibility_strains_plotly}}],

	  		name: 'Strain2', 
	  		marker: {
	  			color: '#0d233a',
	  			line:{
					width: 1,
					color:'black',
				},
			}, 
	  		type: 'bar',
	  		orientation: 'h'
	  	};
	  
	  	var data = [trace1, trace2];
	  
	  	var layout = {
	  		width:  600,
			height: 900,
			
			shapes: [{
				type: 'line',
    			x0: '100',
    			y0: 0,
				x1: '100',
				yref: 'paper',
				y1: 1,
				line: {
					color: 'grey',
					width: 1.5,
					dash: 'dot'
				}
			}],
	  		xaxis: {
	  			automargin: true,
				tickfont: {
	  				size: 14, 
	  			},
				range: [0, 100],
				title:'Strain Compatibility [in %]',
				titlefont: {
	  				size: 18, 
	  			},
			}, 
	  		yaxis: {
				automargin: true,
	  			autorange: 'reversed',
				titlefont: {
	  				size: 18, 
	  			},
	  			tickfont: {
	  				<!-- family: 'Courier New, monospace', -->	
	  				size: 14, 
		  		},
				ticks: 'outside',
				tickcolor: 'rgba(0,0,0,0)',
	  		}, 
			margin: {
			    l: 0,
			    r: 20,
			    b: 100,
			    t: 0,
			    pad: 0,
			},
	  		showlegend: true,
	  		legend: {
	  			"orientation": "h",
	  			x: 0.2, 
	  			y: 1.05, 
				traceorder: 'normal',
	  			font: {
        			size: 16,
        			color: 'black',
      			},	
			},
			
	   		barmode: 'stack', 
	  		bargap: 0.15, 
	  		bargroupgap: 0.1
	  	};
	  
		var optionsHybrids= {
                        displaylogo: false,
						
						modeBarButtonsToRemove:
							['zoom2d',
							'zoom',
							'pan',
							'sendDataToCloud',
							'pan2d',
							'boxSelect',
							'resetScale2d',
							'hoverClosestCartesian',
							'hoverCompareCartesian',
							'toggleSpikelines']
                        ,
                        toImageButtonOptions: {
							filename: 'reStrainingOrder Hybrid Summary',
							width: 600,  <!-- width: stacksDivPercentage._fullLayout.width, height: stacksDivPercentage._fullLayout.height -->
							height: 900,
							format: 'png'
                        }
                };

	  	Plotly.newPlot('Hybrid_Allelic_Ratio_Plot', data, layout, optionsHybrids);

	
	</script>

	<footer>
		<a style="float:right;" href="https://www.bioinformatics.babraham.ac.uk/">
		  {{bioinf_logo_goes_here}}
		  This will be replaced with the logo of Babraham Bioinf
		  {{bioinf_logo_goes_here}}
		</a>

		<p>Analysis produced by <a href="https://github.com/FelixKrueger/reStrainingOrder"><strong>reStrainingOrder</strong></a> (version {{report_version}}) - a Tool to assist with Mouse Strain Identification</p>
		<p>Report graphs rendered using <a href="https://plot.ly/">plot.ly</a>, (v1.48.3), design last changed 21 July 2019</p>
	</footer>

	</div>
</html>
