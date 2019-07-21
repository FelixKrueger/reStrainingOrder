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
                  }
                  .container {
                          margin:0 auto;
                          max-width:1200px;
                  }
                  .header h1,
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
						width:600px;
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
				.error-msg h3 { margin: 0; }
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

	<h2>Strain Compatibility Confidence Level</h2>
	<table>
		<tbody>
			<tr>
				<td style="vertical-align:top">
					<table class="data">
						{{strain_compatibility_confidence}}
					</table>			
				</td>
				<td>
					<div id="Confidence_Plot" class="bargraph"> <!--Plotly chart will be drawn inside this DIV --> </div>
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
					<div id="Strand_Plot" class="bargraph"> <!--Plotly chart will be drawn inside this DIV --> </div>
				</td>
			</tr>
		</tbody>
	</table>

	<hr>

	<!-- This section is optional -->
	<h2>Potential Hybrid Strains - Allelic Ratios</h2>
	<table>
		<tbody>
			<tr>
				<td style="vertical-align:top">
					<table class="data">
						{{allelic_ratio_content}}
					</table>			
				</td>
				<td>
					<div id="Strand_Plot" class="bargraph"> <!--Plotly chart will be drawn inside this DIV --> </div>
				</td>
			</tr>
		</tbody>
	</table>


	<!-- This section is optional -->

	
	<!-- ######################### PLOT.LY plotting code below ################################################################# -->

	<!-- Alignment Stats Plot -->	
	<script>
		<!-- JAVASCRIPT CODE GOES HERE -->
		
		var data = [{
			values: [{{alignment_stats_plotly}}],
			labels: ['Unique Alignments', 'No Alignment', 'Multiple Alignments', 'No Genomic Sequence'],
			type: 'pie',
			name: 'Alignment Statistics',
			hoverinfo: 'label+value+percent',
			direction: 'clockwise',
			pull: [0.05,0,0,0],
			sort: false,
			marker: {
    			colors: ['#0d233a', '#2f7ed8','#8bbc21','#910000','#1aadce','#492970','#f28f43','#77a1e5','#c42525','#a6c96a'],
    			line:{
					width: 1,
					color:'black',
				},
  			},
		}];

		var layout = {
			margin: {
			    l: 50,
			    r: 0,
			    b: 0,
			    t: 50,
			    pad: 0,
			},
			font: {
        		size: 16,
      		},
			<!-- paper_bgcolor: '#7f7f7f', --> 
  			<!-- plot_bgcolor: '#c7c7c7',  -->
			height: 450,
			width:  600,
			showlegend: true,
  			legend: {
    			x: 0.9,
    			y: 1.05,
  				font: {
        			size: 14,
        			color: 'black',
      			},	
  			},
		};

		Plotly.newPlot('Strain_Scores_Report', data, layout, {displaylogo: false}, {modeBarButtonsToRemove: ['toImage',
					'sendDataToCloud',
					'resetScale2d',
					'hoverClosestCartesian',
                    'hoverCompareCartesian',
                    'toggleSpikelines']},
					);
	</script>

	<!-- Duplication Plot (Donut Plot)-->	
	<script>
		
		var data = [{
			values: [{{duplication_stats_plotly}}],
			labels: ['Unique Alignments', 'Duplicate Alignments' ],
			textinfo: ['Unique Alignments', 'Duplicate Alignments' ],
			name: 'Degree of Duplication',
			hoverinfo: 'label+value+percent',
			hole: .3,	
			type: 'pie',
			direction: 'clockwise',
			pull: [0.05,0],
			sort: false,
			marker: {
				line:{
					width: 1,
					color:'black',
				},
				colors: ['#0d233a', '#2f7ed8'],
			}
	 	}];

		var layout = {
			<!-- paper_bgcolor: '#7f7f7f', -->
  			<!-- plot_bgcolor: '#c7c7c7',  -->
  			font: {
        		size: 16,
        		color: 'white',
      		},
      		margin: {
			    l: 0,
			    r: 0,
			    b: 0,
			    t: 0,
			    pad: 0,
			},
			showlegend: true,
  			legend: {
    			x: .25,
    			y: -0.05,
  			        "orientation": "h",
  				font: {
        			size: 14,
        			color: 'black',
      			},	
  			},
		};

		Plotly.newPlot('Duplication_Plot', data, layout, {displaylogo: false}, {modeBarButtonsToRemove: ['toImage',
					'sendDataToCloud',
					'resetScale2d',
					'hoverClosestCartesian',
                    'hoverCompareCartesian',
                    'toggleSpikelines']});
	</script>

	<!-- Cytosine Methylation Plot -->
	<script>
		var data = [
		  {
			x: ['CpG context', 'CHG context', 'CHH context'],
			y: [{{cytosine_methylation_plotly}}],
			type: 'bar',
			marker: {
				color: ['#0d233a', '#2f7ed8','#2f7ed8'],
	 			<!-- colors: ['#0d233a', '#2f7ed8','#8bbc21','#910000','#1aadce','#492970','#f28f43','#77a1e5','#c42525','#a6c96a'], -->
   				line: {
     				color: 'black',
      				width: 1
    			},
   			}	
		  }
		];
		
		var layout = {
			<!-- paper_bgcolor: '#7f7f7f', -->
  			<!-- plot_bgcolor: '#c7c7c7',  -->
  			width: 600,
  			margin: {
			    l: 100,
			    r: 50,
			    b: 50,
			    t: 50,
			    pad: 0,
			},
			yaxis: {
				range: [0, 100],
				title:' % Methylation',
			}
		};

		Plotly.newPlot('Cytosine_Methylation', data, layout, {displaylogo: false}, {modeBarButtonsToRemove: ['toImage',
					'sendDataToCloud',
					'resetScale2d',
					'hoverClosestCartesian',
                    'hoverCompareCartesian',
                    'toggleSpikelines']});
	</script>

	<!-- Cytosine Methylation Plot post duplication-->
	<script>

		var data = [
		  	{
				x: ['CpG context', 'CHH context', 'CHG context'],
				y: [{{cytosine_methylation_post_duplication_plotly}}],
				type: 'bar',
				marker: {
					color: ['#0d233a', '#2f7ed8','#2f7ed8'],
	 				line: {
     					color: 'black',
      					width: 1
    				},
				}	
		  	}
		];

		var layout = {
			width: 600,
			<!-- paper_bgcolor: '#7f7f7f', -->
  			<!-- plot_bgcolor: '#c7c7c7',  -->
  			margin: {
			    l: 100,
			    r: 50,
			    b: 50,
			    t: 50,
			    pad: 0,
			},
			yaxis: {
				range: [0, 100],
				title:' % Methylation',
			}
		};

		Plotly.newPlot('Cytosine_Methylation_postDuplication', data, layout, 
					{displaylogo: false}, 
					{modeBarButtonsToRemove: 
					['toImage',
					'sendDataToCloud',
					'resetScale2d',
					'hoverClosestCartesian',
                    'hoverCompareCartesian',
                    'toggleSpikelines']
					}
		);
	</script>

	<!-- Strand Alignment Plot-->
	<script>
		var data = [
		  {
			x: ['OT', 'CTOT', 'CTOB', 'OB'],
			<!-- y: [49, 1, 2, 48], hardcoded for testing purposes -->
			y: [{{strand_alignment_plotly}}],
			type: 'bar',
			marker: {
				color: ['#0d233a', '#2f7ed8','#2f7ed8','#0d233a'],
	 			line: {
     				color: 'black',
      				width: 1
    			},
   			}	
		  }
		];

		var layout = {
			<!-- paper_bgcolor: '#7f7f7f',  -->
  			<!-- plot_bgcolor: '#c7c7c7',   -->
  			margin: {
			    l: 100,
			    r: 50,
			    b: 50,
			    t: 15,
			    pad: 5,
			},
			yaxis: {
				title: 'Number of Alignments',
				zeroline:true, 
				<!-- hoverformat: '.2r', -->
			}
		};


		Plotly.newPlot('Strand_Plot', data, layout, {displaylogo: false}, {modeBarButtonsToRemove: ['toImage',
					'sendDataToCloud',
					'resetScale2d',
					'hoverClosestCartesian',
                    'hoverCompareCartesian',
                    'toggleSpikelines']});
	</script>

	<!-- Nucleotide Stats Plot-->
	<script>
	  	var trace1 = {
	  		x: [{{nucleo_sample_x}}], 
	  		y: [{{nucleo_sample_y}}],
	  
	  		name: 'Percent Sample', 
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
	  		<!-- x: [20, 14, 23, 7], testing -->
	  		<!-- y: ['A', 'T', 'C', 'G'], testing -->
	 		x: [{{nucleo_genomic_x}}], 
	  		y: [{{nucleo_genomic_y}}],

	  		name: 'Percent Genomic', 
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
	  		height: 700,
	  		width:  700,
 
	  		xaxis: {
	  			tickfont: {
	  				size: 14, 
	  			}
	  		}, 
	  		yaxis: {
	  			title: '(Di-)Nucleotide',
	  			autorange: 'reversed',
	  			titlefont: {
	  				size: 18, 
	  			},
	  			tickfont: {
	  				<!-- family: 'Courier New, monospace', -->	
	  				size: 14, 
		  		}
	  		}, 
	  		showlegend: true,
	  		legend: {
	  			"orientation": "h",
	  			x: 0.15, 
	  			y: -0.05, 
	  			font: {
        			size: 16,
        			color: 'black',
      			},	
			},
	  
	  		barmode: 'group', 
	  		bargap: 0.15, 
	  		bargroupgap: 0.1
	  	};
	  
	  	Plotly.newPlot('nucleo_plot', data, layout, {displaylogo: false},{modeBarButtonsToRemove: ['toImage',
					'sendDataToCloud',
					'resetScale2d',
					'hoverClosestCartesian',
                    'hoverCompareCartesian',
                    'zoom',
                    'toggleSpikelines']});
	
	</script>

	<!-- M-bias Plot 1 -->
	<script>
		var trace1 = {
	  		<!-- x: [1, 2, 3,4,5,6,7,8,9,10], just for testing -->
	  		<!-- y: [40, 50, 60,50,35,40,45,61,55,33], just for testing -->
	  		<!-- colors: ['#0d233a', '#2f7ed8','#8bbc21','#910000','#1aadce','#492970','#f28f43','#77a1e5','#c42525','#a6c96a'], -->
	  		<!-- colors: [ '#CCF0E1','#EDD3A8','#69798A','#21BCA2','#F29D13','#0d233a','#f28f43','#77a1e5','#c42525','#a6c96a'], -->
	  		x: [{{mbias1_CpG_meth_x}}],
	  		y: [{{mbias1_CpG_meth_y}}],
	  		name: 'CpG methylation',
	  		line: {
	  			width: 5,
	  			color: '#0d233a',
	  		},
	  		type: 'scatter',
	  		yaxis: 'y',
		};

		var trace2 = {
			x: [{{mbias1_CHG_meth_x}}],
	  		y: [{{mbias1_CHG_meth_y}}],
	  		name: 'CHG methylation',
	  		yaxis: 'y',
	  		type: 'scatter',
	  		line: {
	  			width: 5,
	  			color: '#F29D13',
	  		},
		};

		var trace3 = {
			x: [{{mbias1_CHH_meth_x}}],
	  		y: [{{mbias1_CHH_meth_y}}],
	  		name: 'CHH methylation',
	  		yaxis: 'y',
	  		type: 'scatter',
	  		line: {
	  			width: 5,
	  			color: '#21BCA2',
	  		},
		};

		var trace4 = {
	  		
	  		x: [{{mbias1_CpG_coverage_x}}],
	  		y: [{{mbias1_CpG_coverage_y}}],
	  		name: 'CpG total calls',
	  		type: 'scatter',
	  		yaxis: 'y2',
	  		opacity: 0.4,
	  		line: {
	  			width: 1.5,
	  			color: '#0d233a',
	  		},
		};

		var trace5 = {
			x: [{{mbias1_CHG_coverage_x}}],
	  		y: [{{mbias1_CHG_coverage_y}}],
	  		name: 'CHG total calls',
	  		yaxis: 'y2',
	  		type: 'scatter',
	  		opacity: 0.4,
	  		line: {
	  			width: 1.5,
	  			color: '#F29D13',
	  		},
		};

		var trace6 = {
			x: [{{mbias1_CHH_coverage_x}}],
	  		y: [{{mbias1_CHH_coverage_y}}],
	  		name: 'CHH total calls',
	  		yaxis: 'y2',
	  		type: 'scatter',
	  		opacity: 0.4,
	  		line: {
	  			width: 1.5,
	  			color: '#21BCA2',
	  		},
		};

		var data = [trace1, trace2, trace3, trace4, trace5, trace6];

		var layout = {
	  		title: 'Read 1',
	  		titlefont: {
			      	size: 24,
			      	color: 'black'
			},
	  		yaxis: {
	  			title: '% Methylation',
	  			range: [0, 100],
	  			visible: true,
	  			titlefont: {
			      	size: 18,
			      	color: 'black'
			    },
			    showline: true,
	  		},
	  		width: 1200,
	  		height: 600,
	  		yaxis2: {
	    		title: '# Methylation Calls',
	    		tickfont: {color: 'black'},
	    		overlaying: 'y',
	    		side: 'right',
	    		separatethousands: 'y',
	    		zeroline: false,
	    		visible: true,
	    		showgrid: false,
	    		showline: true,
	    		titlefont: {
			      	size: 18,
			      	color: 'black'
			    },
	  		},
	  		xaxis: {
	  			title: 'Position in Read [bp]',
	  			titlefont: {
			      	size: 18,
			      	color: 'black'
			    },
			    showline: true,
	  		},
	  		showlegend: true,
		  	legend: {
		    	x: 0.85,
		    	y: 1.2
		  	},
		};

		var options = {
			displaylogo: false,
			modeBarButtonsToRemove:
				['zoom2d',
				'sendDataToCloud',				
				'pan', 
				'pan2d',
				'resetScale2d',
				'hoverClosestCartesian',
				'hoverCompareCartesian',
				'toggleSpikelines']
			,
			toImageButtonOptions: {
				filename: 'Bismark M-bias Read 1',
				width: 1600,
				height: 600,
				format: 'png'
			}
		};
		
		Plotly.newPlot('mbias1_plot', data, layout, options);
	</script>

	<!-- M-bias Plot 2-->
	<script>
		var trace1 = {
	  		<!-- x: [1, 2, 3,4,5,6,7,8,9,10], just for testing -->
	  		<!-- y: [40, 50, 60,50,35,40,45,61,55,33], just for testing -->

	  		x: [{{mbias2_CpG_meth_x}}],
	  		y: [{{mbias2_CpG_meth_y}}],
	  		name: 'CpG methylation',
	  		type: 'scatter',
	  		yaxis: 'y',
	  		line: {
	  			width: 5,
	  			color: '#0d233a',
	  		},
		};

		var trace2 = {
			x: [{{mbias2_CHG_meth_x}}],
	  		y: [{{mbias2_CHG_meth_y}}],
	  		name: 'CHG methylation',
	  		yaxis: 'y',
	  		type: 'scatter',
	  		line: {
	  			width: 5,
	  			color: '#F29D13',
	  		},
		};

		var trace3 = {
			x: [{{mbias2_CHH_meth_x}}],
	  		y: [{{mbias2_CHH_meth_y}}],
	  		name: 'CHH methylation',
	  		yaxis: 'y',
	  		type: 'scatter',
	  		line: {
	  			width: 5,
	  			color: '#21BCA2',
	  		},
		};

		var trace4 = {		
	  		x: [{{mbias2_CpG_coverage_x}}],
	  		y: [{{mbias2_CpG_coverage_y}}],
	  		name: 'CpG total calls',
	  		type: 'scatter',
	  		yaxis: 'y2',
	  		opacity: 0.4,
	  		line: {
	  			width: 1.5,
	  			color: '#0d233a',
	  		},
		};

		var trace5 = {
			x: [{{mbias2_CHG_coverage_x}}],
	  		y: [{{mbias2_CHG_coverage_y}}],
	  		name: 'CHG total calls',
	  		yaxis: 'y2',
	  		type: 'scatter',
	  		opacity: 0.4,
	  		line: {
	  			width: 1.5,
	  			color: '#F29D13',
	  		},
		};

		var trace6 = {
			x: [{{mbias2_CHH_coverage_x}}],
	  		y: [{{mbias2_CHH_coverage_y}}],
	  		name: 'CHH total calls',
	  		yaxis: 'y2',
	  		type: 'scatter',
	  		opacity: 0.4,
	  		line: {
	  			width: 1.5,
	  			color: '#21BCA2',
	  		},
		};

		var data = [trace1, trace2, trace3, trace4, trace5, trace6];

		var layout = {
	  		title: 'Read 2',
	  		titlefont: {
			      	size: 24,
			      	color: 'black'
			},
	  		yaxis: {
	  			title: '% Methylation',
	  			range: [0, 100],
	  			titlefont: {
			      	size: 18,
			      	color: 'black'
			    },
			    showline: true,
	  		},
	  		width: 1200,
	  		height: 600,
	  		yaxis2: {
	    		title: '# Methylation Calls',
	    		overlaying: 'y',
	    		side: 'right',
	    		separatethousands: 'y',
	    		zeroline: false,
	    		showgrid: false,
	    		showline: true,
	    		titlefont: {
			      	size: 18,
			      	color: 'black'
			    },
	  		},
	  		xaxis: {
	  			title: 'Position in Read [bp]',
	  			showline: true,	
	  			titlefont: {
			      size: 18,
			      color: 'black'
			    },
	  		},
	  		showlegend: true,
		  	legend: {
		    	x: 0.85,
		    	y: 1.2,
		  	},
		};

		
		var options2 = {
			displaylogo: false,
			modeBarButtonsToRemove:
				['zoom2d',
				'sendDataToCloud',				
				'pan', 
				'pan2d',
				'resetScale2d',
				'hoverClosestCartesian',
				'hoverCompareCartesian',
				'toggleSpikelines']
			,
			toImageButtonOptions: {
				filename: 'Bismark M-bias Read 2',
				width: 1600,
				height: 600,
				format: 'png'
			}
		};
		
		Plotly.newPlot('mbias2_plot', data, layout, options2);

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
