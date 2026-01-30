from cwest_polymer import PolyGraph
from cwest_polymer import FractionalMRDataset
from cwest_polymer import fmr_transforms as transforms
from piblin.data import Measurement, MeasurementSet
import os


def test_fractional_mr_networkgraph(mass, rt, abundance, repeat_units,test_dir):
    fmr_graph = PolyGraph()

    measurementset = MeasurementSet(
        measurements=[
            Measurement(
                datasets=[
                    FractionalMRDataset(
                        mass=mass,
                        rt=rt,
                        abundance=abundance
                    )
                ],
                conditions={'file_name': 'test_file'}
            )
        ],
        merge_redundant=False
    )

    fmr_transform = transforms.FractionalMRTransform.create(
        repeat_units=repeat_units,
        default_list=False,
        fractional_values=1
    )

    fmr_cluster = transforms.ClusterTransform.create(
        mz_tol=0.001,
        ppm_tol=10,
        min_samples=3
    )

    pipeline = fmr_transform + fmr_cluster

    result = pipeline(measurementset)

    fmr_graph.add_measurements(result.measurements)
    assert fmr_graph.graph.number_of_nodes() == len(mass)
    assert fmr_graph.graph.number_of_edges() == len(mass) * 2

    fig = fmr_graph.plot_graph_with_plotly()
    assert len(fig.data) == 3
    assert len(fig.data[1]['text']) == len(mass) * 2

    fig1 = fmr_graph.plot_spectrum_with_plotly()
    assert len(fig1.data) == 5
    assert len(fig.data[1]['text']) == len(mass) * 2

    temp_dir = os.path.join(test_dir, "temp_fmr_graph.html")
    try:
        fmr_graph.save_plot_to_html(file_path=temp_dir)
        with open(temp_dir, 'r', encoding='utf-8') as file:
            temp_content = file.read()
        os.remove(temp_dir)
    except Exception as e:
        if os.path.isfile(temp_dir):
            os.remove(temp_dir)
        raise e
    assert temp_content.startswith('<html>\n<head>')
    assert temp_content.endswith('</body>\n</html>')
