package krigingsTests;

import java.util.HashMap;
import java.util.Iterator;
import java.util.Set;

import org.geotools.coverage.grid.GridCoverage2D;
import org.geotools.coverage.grid.GridGeometry2D;
import org.geotools.data.simple.SimpleFeatureCollection;
import org.geotools.filter.text.cql2.CQL;
import org.geotools.resources.coverage.CoverageUtilities;
import org.jgrasstools.gears.io.rasterreader.OmsRasterReader;
import org.jgrasstools.gears.io.rasterwriter.OmsRasterWriter;
import org.jgrasstools.gears.io.shapefile.OmsShapefileFeatureReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorReader;
import org.jgrasstools.gears.io.timedependent.OmsTimeSeriesIteratorWriter;
import org.opengis.feature.simple.SimpleFeature;
import org.opengis.filter.Filter;
import java.awt.geom.Point2D;

import com.vividsolutions.jts.geom.Coordinate;
import com.vividsolutions.jts.geom.Geometry;

import org.jgrasstools.hortonmachine.utils.HMTestCase;


import krigingsRasterCase.Krigings;


public class TestKrigingRasterCase extends HMTestCase {


	/**
	 * Run the kriging models.
	 *
	 * <p>
	 * This is the case which all the station have the same value.
	 * </p>
	 * @throws Exception
	 * @throws Exception
	 */

	public void testKriging() throws Exception {
		//
		String stationIdField = "ID_P";

		OmsShapefileFeatureReader stationsReader = new OmsShapefileFeatureReader();
		stationsReader.file = "resources/Input/krigings/RasterCase/final_station.shp";
		stationsReader.readFeatureCollection();
		SimpleFeatureCollection stationsFC = stationsReader.geodata;

		// OmsShapefileFeatureReader interpolatedPointsReader = new OmsShapefileFeatureReader();
		// interpolatedPointsReader.file = puntiFile.getAbsolutePath();
		// interpolatedPointsReader.readFeatureCollection();

		OmsTimeSeriesIteratorReader reader = new OmsTimeSeriesIteratorReader();
		reader.file = "resources/Input/krigings/RasterCase/P1994_2013.csv";
		reader.idfield = "ID_P";
		reader.tStart = "1994-01-01 00:00";
		reader.tTimestep = 60;
		 reader.tEnd = "1994-01-01 00:00";
		reader.fileNovalue = "-999";
		reader.initProcess();

		Krigings kriging = new Krigings();
		kriging.pm = pm;

		
		OmsRasterReader demReader = new OmsRasterReader();
		demReader.file = "resources/Input/krigings/RasterCase/dem.asc";
		demReader.fileNovalue = -9999.0;
		demReader.geodataNovalue = Double.NaN;
		demReader.process();
		GridCoverage2D dem = demReader.outRaster;

		
		
		kriging.inGridCoverage2D = dem;

		kriging.inStations = stationsFC;
		kriging.fStationsid = stationIdField;

		kriging.range = 123537.0;
		kriging.nugget = 0.0;
		kriging.sill= 1.678383;
		kriging.pSemivariogramType="linear";


		while( reader.doProcess ) {
			reader.nextRecord();
			HashMap<Integer, double[]> id2ValueMap = reader.outData;
			kriging.inData = id2ValueMap;
			kriging.executeKriging();
			/*
			 * Extract the result.
			 */



			GridCoverage2D krigingRaster = kriging.outGrid;
			
			OmsRasterWriter writerRaster = new OmsRasterWriter();
			writerRaster.inRaster = krigingRaster;
			writerRaster.file = "resources/Output/krigings/RasterCase/krigings.asc";
			writerRaster.process();

			
		}

		reader.close();
		// writer.close();
	}
}
