<HTML>
<HEAD>

	<LINK rel="Stylesheet" href="main.css" type="text/css" media="screen">	
</HEAD>


<?php

try{
$con = new PDO("sqlite:../analysis_manager.db3");

$sql = "SELECT filename FROM features_analysis_plots WHERE format_id = 'output_web_raster';";

$res = $con->query($sql);

$plot_id = 1;
foreach($res as $row){
	echo "<DIV class=\"plotdiv\" align=\"center\">\n";
	echo "<IMG src=\"../", $row['filename'], "\" class=\"plot\" id=\"plot_", $plot_id, "\"/>\n";
	echo "</DIV>\n";
	echo "<DIV class=\"plot_separator\"></DIV>\n\n";
	
	$plot_id += 1;
}

$con = NULL;

} catch (PDOException $e) {
	print 'Exception : '.$e->getMessage();
}
?>
</HTML>
