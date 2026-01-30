"""
This module contains the PolyGraph class, which is used to create and visualize polymer graphs
"""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots

import networkx as nx
from .fmr_datasets import FractionalMRDataset as FMR
from .. import fmr_parameters as p
from piblin.data import Measurement
from typing import Tuple, List
from pathlib import Path


class PolyGraph:
    """Class for creating and visualizing polymer graphs based on polymeric analysis results.

    This class provides methods to construct polymer graphs from datasets, add measurements, and visualize them using
    interactive plots. It also includes properties to access the graph and its nodes.

    Methods
    -------
    add_fmr_results(fmr_result: FMR, ru_label: Tuple)
        Adds fractional MR results to the polymer graph.
    add_measurements(measurements: List[Measurement])
        Adds multiple measurements to the polymer graph.
    plot_graph_with_plotly(node_size=10, node_color='blue', edge_label_size=10)
        Plots the polymer graph using Plotly.
    plot_spectrum_with_plotly(node_size: float = 10, node_color: str = 'blue')
        Plots the spectrum view of the polymer graph using Plotly.
    save_plot_to_html(file_path: Path)
        Saves the polymer graph plot to an HTML file with interactive features.

    Attributes
    ----------
    graph : nx.Graph
        The NetworkX graph representing the polymer structure.
    """

    def __init__(self):
        self.graph = nx.Graph()

    def add_fmr_results(self, fmr_result: FMR, ru_label: Tuple):
        ru_label = f"{ru_label[0]}_{ru_label[1]:.3f}"
        data_array = fmr_result.data_arrays
        data_array_names = fmr_result.data_array_names
        cmpd_list = data_array[data_array_names.index(p.CMPD_LABEL)]
        mass_list = data_array[data_array_names.index(p.MASS_LABEL)]
        abundance = data_array[data_array_names.index(p.ABUNDANCE_LABEL)]
        clusters = data_array[data_array_names.index(p.CLUSTER_LABEL)]
        # nodes
        for n, mass, abund in zip(cmpd_list.tolist(), mass_list.tolist(), abundance.tolist()):
            if n in self.graph.nodes:
                continue
            self.graph.add_node(n, pos=(mass, abund))
        # edges
        unique_clusters = [x for x in set(clusters) if x not in p.UNCLUSTERED_LABELS]
        for cluster in unique_clusters:
            cluster_numbers = cmpd_list[np.where(clusters == cluster)[0]].tolist()
            cluster_masses = mass_list[np.where(clusters == cluster)[0]].tolist()
            cluster_masses.sort()
            for n, g1 in zip(cluster_numbers, cluster_masses):
                for m, g2 in zip(cluster_numbers, cluster_masses):
                    if g1 != g2 and not self.graph.has_edge(n, m):
                        label = ru_label
                        self.graph.add_edge(n, m, label=label)

    def add_measurements(self, measurements: List[Measurement]):
        for measurement in measurements:
            dataset = measurement.datasets[0]
            if (dataset.number_of_points() == 0) or (isinstance(dataset, FMR) is False) or (
                    measurement.num_datasets != 1):
                print(f"Measurement does not a valid FMR dataset.\n{measurement.conditions.get('file_name', '')}")
                continue
            ru_info = measurement.details['repeat_unit_information']
            if isinstance(measurement.datasets[0], FMR):
                self.add_fmr_results(measurement.datasets[0], ru_info)

    def plot_graph_with_plotly(self, node_size=10, node_color='blue', edge_label_size=10):
        # Generate positions for nodes
        pos = nx.spring_layout(self.graph, seed=42)  # Fixed seed for reproducibility
        nx.set_node_attributes(self.graph, pos, 'spring_pos')

        # Create edge traces
        edge_x = []
        edge_y = []
        edge_labels = []
        edge_label_positions = []
        for n0, n1, details in self.graph.edges(data=True):
            x0, y0 = pos[n0]
            x1, y1 = pos[n1]
            label = details.get('label', '')
            edge_x.extend([x0, x1, None])
            edge_y.extend([y0, y1, None])
            edge_labels.append(label)
            # center position for edge label
            edge_label_positions.append(((x0 + x1) / 2, (y0 + y1) / 2))

        edge_trace = go.Scatter(
            x=edge_x, y=edge_y,
            line=dict(width=0.5, color='#888'),
            hoverinfo='none',
            mode='lines'
        )

        # Create edge label traces
        if len(edge_labels) > 0:
            edge_label_x, edge_label_y = zip(*edge_label_positions)
            edge_label_trace = go.Scatter(
                x=edge_label_x, y=edge_label_y,
                mode='text',
                text=edge_labels,
                textfont=dict(size=edge_label_size, color='red'),
                hoverinfo='none'
            )

        # Create node traces
        node_x = []
        node_y = []
        node_text = []
        for node, details in self.graph.nodes(data=True):
            mass, abund = details['pos']
            spring_x, spring_y = details['spring_pos']
            label = f"{mass:.3f}"
            node_x.append(spring_x)
            node_y.append(spring_y)
            node_text.append(label)  # Add node labels as hover text

        node_trace = go.Scatter(
            x=node_x, y=node_y,
            mode='markers+text',  # Add text mode to display labels
            hoverinfo='text',
            text=node_text,  # Display node labels
            textposition='top center',  # Position labels above the nodes
            marker=dict(
                size=node_size,
                color=node_color,
                line_width=2
            )
        )

        # Combine traces into a figure
        if len(edge_labels) > 0:
            traces = [edge_trace, edge_label_trace, node_trace]
        else:
            traces = [edge_trace, node_trace]
        fig = go.Figure(data=traces,
                        layout=go.Layout(
                            showlegend=False,
                            hovermode='closest',
                            margin=dict(b=0, l=0, r=0, t=0),
                            xaxis=dict(showgrid=False, zeroline=False),
                            yaxis=dict(showgrid=False, zeroline=False)
                        ))

        return fig

    def plot_spectrum_with_plotly(self, node_size: float = 10, node_color: str = 'blue'):
        # Generate positions for nodes using spring layout
        spring_pos = nx.spring_layout(self.graph, seed=42)  # Fixed seed for reproducibility
        nx.set_node_attributes(self.graph, spring_pos, 'spring_pos')

        # Create traces for spring layout graph
        spring_edge_x, spring_edge_y = [], []
        spring_edge_labels = []
        for n1, n0, details in self.graph.edges(data=True):
            x0, y0 = spring_pos[n0]
            x1, y1 = spring_pos[n1]
            spring_edge_x.extend([x0, x1, None])
            spring_edge_y.extend([y0, y1, None])
            spring_edge_labels.extend([details.get('label', ''), '', None])

        # Create traces for nodes graphs
        spring_node_x, spring_node_y, spring_node_text = [], [], []
        stored_node_x, stored_node_y, stored_node_text = [], [], []
        for n0, details in self.graph.nodes(data=True):
            mass, abund = details['pos']
            spring_x, spring_y = details['spring_pos']
            label = f"{mass:.3f}"

            spring_node_x.append(spring_x)
            spring_node_y.append(spring_y)
            spring_node_text.append(label)

            stored_node_x.append(mass)
            stored_node_y.append(abund)
            stored_node_text.append(label)

        # create scatter plots
        spring_edge_trace = go.Scatter(
            x=spring_edge_x, y=spring_edge_y,
            line=dict(width=0.5, color='#888'),
            hoverinfo='none',
            mode='lines'
        )

        # create separate trace for edge labels
        spring_edge_label_trace = go.Scatter(
            x=[(x0 + x1) / 2 for x0, x1, _ in zip(spring_edge_x[::3], spring_edge_x[1::3], spring_edge_x[2::3])],
            y=[(y0 + y1) / 2 for y0, y1, _ in zip(spring_edge_y[::3], spring_edge_y[1::3], spring_edge_y[2::3])],
            mode='text',
            text=[label for label in spring_edge_labels if label is not None],
            textfont=dict(size=8, color='rgba(0, 0, 255, 0.5)'),
            hoverinfo='none'
        )

        spring_node_trace = go.Scatter(
            x=spring_node_x, y=spring_node_y,
            mode='markers',
            hoverinfo='text',
            text=spring_node_text,
            marker=dict(
                size=node_size,
                color=node_color,
                line_width=2
            ),
            name='Spring Layout'
        )

        stored_lines_trace = go.Scatter(
            x=[val for x in stored_node_x for val in [x, x, None]],  # Repeat each x value twice
            y=[val for y in stored_node_y for val in [y, 0, None]],  # Connect each point to y=0
            mode='lines',
            line=dict(color=node_color, width=0.5),
            hoverinfo='none',
            showlegend=False,
            name='Stored Lines'
        )

        stored_positions_trace = go.Scatter(
            x=stored_node_x,
            y=stored_node_y,
            mode='markers',
            hoverinfo='text',
            text=stored_node_text,
            marker=dict(
                size=node_size,
                color=node_color,
                line_width=2
            ),
            name='Stored Positions'
        )

        # Create subplots for both layouts
        fig = make_subplots(rows=1, cols=2, subplot_titles=("Network Graph", "Spectrum View"))

        # Add traces to subplots
        fig.add_trace(spring_edge_trace, row=1, col=1)
        fig.add_trace(spring_node_trace, row=1, col=1)

        fig.add_trace(stored_lines_trace, row=1, col=2)
        fig.add_trace(stored_positions_trace, row=1, col=2)

        fig.add_trace(spring_edge_label_trace, row=1, col=1)

        # Update layout for interactivity
        fig.update_layout(
            title='Polymer Spectrum with Network Overlay',
            xaxis2_title='Mass (m/z)',
            yaxis2_title='Abundance',
            dragmode='select',
            modebar_add=['scrollZoom']
        )

        return fig

    def save_plot_to_html(self, file_path: Path):
        fig = self.plot_spectrum_with_plotly(node_size=4)
        fig.write_html(file_path)

        with open(file_path, 'rb') as f:
            html_content = f.read().decode('utf-8')

        # Insert the linking JavaScript just before the closing </body> tag
        js_code = """
        <script>
        (function() {
            let retries = 0;
            let defNodeSize;
            document.addEventListener('DOMContentLoaded', () => {
                const mainPlot = document.querySelector('.plotly-graph-div');
                if (mainPlot && mainPlot.data && mainPlot.data.length >= 5) {
                    defNodeSize = mainPlot.data[1].marker.size;  // Get initial size from spring plot nodes
                } else {
                    defNodeSize = 4;  // Fallback default size
                }
            });
                function linkPlots() {
                    // Wait for plots and Plotly to be fully loaded
                    const maxRetries = 50; // 5 seconds maximum wait
                    if (typeof Plotly === 'undefined') {
                        console.log(`Waiting for Plotly to load... (attempt ${retries + 1}/${maxRetries})`);
                        if (retries < maxRetries) {
                            retries++;
                            setTimeout(linkPlots, 100);
                        } else {
                            console.error('Plotly failed to load after 5 seconds');
                        }
                        return;
                    }
    
                    var mainPlot = document.querySelector('.plotly-graph-div');
                    if (!mainPlot || !mainPlot.data || mainPlot.data.length < 5) {
                        console.log('Main plot found:', !!mainPlot);
                        console.log('Plot data ready:', mainPlot ? mainPlot.data.length : 0);
                        if (retries < maxRetries) {
                            retries++;
                            setTimeout(linkPlots, 100);
                        } else {
                            console.error('Fully rendered plot div not found after 5 seconds');
                        }
                        return;
                    }
            
                    var springGraphDiv = mainPlot;  // First subplot is in the main div
                    var storedGraphDiv = mainPlot;  // Second subplot is in the main div
                    
                    console.log('Found plot divs:', springGraphDiv.id, storedGraphDiv.id);
            
                    function updateColors(selectedMasses) {
                        console.log('Updating colors for masses:', selectedMasses);
                        try {
                            // Update spring layout plot (traces 0 and 1 in first subplot)
                            var springColors = mainPlot.data[1].text.map(mass => 
                                selectedMasses.has(mass) ? 'red' : 'blue'
                            );
                            var springSizes = mainPlot.data[1].text.map(mass => 
                                selectedMasses.has(mass) ? defNodeSize + 5 : defNodeSize
                            );
                            var storedColors = mainPlot.data[3].text.map(mass => 
                                selectedMasses.has(mass) ? 'red' : 'blue'
                            );
                            var storedSizes = mainPlot.data[3].text.map(mass => 
                                selectedMasses.has(mass) ? defNodeSize + 5 : defNodeSize
                            );
                    
                            // Update spring plot nodes (trace index 1)
                            Plotly.restyle(mainPlot, {
                                'marker.color': [springColors],
                                'marker.size': [springSizes]
                            }, [1]);
                    
                            // Update stored plot nodes (trace index 3)
                            Plotly.restyle(mainPlot, {
                                'marker.color': [storedColors],
                                'marker.size': [storedSizes]
                            }, [3]);
                    
                        } catch (error) {
                            console.error('Error updating colors:', error);
                            console.error('Error details:', error.stack);
                        }
                    }
            
                    function resetColors() {
                        console.log('Resetting colors and sizes');
                        try {
                            Plotly.restyle(springGraphDiv, {
                                'marker.color': ['blue'],
                                'marker.size': [defNodeSize]
                            }, [1]);
                            Plotly.restyle(storedGraphDiv, {
                                'marker.color': ['blue'],
                                'marker.size': [defNodeSize]
                            }, [3]);
                        } catch (error) {
                            console.error('Error resetting colors and sizes:', error);
                        }
                    }

            
                    springGraphDiv.on('plotly_selected', function(eventData) {
                        console.log('Spring plot selection:', eventData);
                        if (!eventData || eventData.points.length === 0) {
                            resetColors();
                            return;
                        }
                        var selectedMasses = new Set(eventData.points.map(p => p.text));
                        updateColors(selectedMasses);
                    });
            
                    storedGraphDiv.on('plotly_selected', function(eventData) {
                        console.log('Stored plot selection:', eventData);
                        if (!eventData || eventData.points.length === 0) {
                            resetColors();
                            return;
                        }
                        var selectedMasses = new Set(eventData.points.map(p => p.text));
                        updateColors(selectedMasses);
                    });
            
                    springGraphDiv.on('plotly_deselect', resetColors);
                    storedGraphDiv.on('plotly_deselect', resetColors);
            
                    console.log('Plot linking initialized');
                }
        
                // Start the initialization process
                document.addEventListener('DOMContentLoaded', function() {
                    console.log('DOM loaded, starting plot linking...');
                    linkPlots();
                });
            
                // Fallback initialization
                if (document.readyState === 'complete') {
                    console.log('Document already complete, starting plot linking...');
                    linkPlots();
                }
            })();
        </script>
        """

        # Insert the JavaScript code before </body>
        html_content = html_content.replace('</body>', f'{js_code}</body>')
        with open(file_path, 'wb') as f:
            f.write(html_content.encode('utf-8'))
