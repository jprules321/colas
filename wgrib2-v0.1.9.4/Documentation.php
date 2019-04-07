<? $root="Resources/Documentation/opengrads/";  include ("doc/opengrads/header.php"); ?>

<h1> Contents </h1>

<ul>

<li> <a href="#cola"> GrADS Documentation from COLA</a>  </li>
<li> <a href="#udxt"> OpenGrADS User Defined Extensions </a>  </li>
<li> <a href="#gagui"> GrADS with Athena Widgets </a>  </li>
<li> <a href="#bundle"> OpenGrADS Bundle Front-ends </a>  </li>
<li> <a href="#scripts"> Scripts </a>  </li>

</ul>

<h1><a name="cola"></a> GrADS Documentation from COLA</h1>

The HTML version of the GrADS documentation distributed by COLA has
become the standard base documentation for GrADS. Check the
<a href="http://grads.iges.org/grads/gaResources/Documentation/index.html"> on-line </a>
version of this documentation for a more up-to-date version.
<ul>

  <li> The <a href="Resources/Documentation/users.html">User's Guide</a></h2> is the
       fundamental document that provides information about how to use
       GrADS. The four main chapters are General Topics, Analysis
       Topics, Display Topics, and the GrADS Scripting Language

  <li> The <a href="Resources/Documentation/tutorial.html">Tutorial</a> will give you a
       feeling for how to use the basic capabilities of GrADS. This
       sample session takes about 30 minutes to run through. It is
       highly recommended for new users.

  <li> The <a href="Resources/Documentation/gadocindex.html">Index</a> provides a quick
       and easy interface for checking the syntax and usage of any
       GrADS command or function. Subject headings from the User's
       Guide are also listed in the Index. </span></p>

</ul>


<h1> <a name="udxt"</a> OpenGrADS User Defined Extensions  </h1>

Here you find documentation for the OpenGrADS <i> User Defined
Commands </i> and <i>User Defined Functins</i> that extends the
GrADS functionality provided by COLA.

   <ul>

   <li><a href="Resources/Documentation/opengrads/doc/udxt/env/index.html">
                                   Env </a> - Environment Variable Manipulation
   </li>

   <li><a href="Resources/Documentation/opengrads/doc/udxt/fish/index.html">
                                   Fish </a> - Streamfunction and Velocity Potential 
   </li>

   <li><a href="Resources/Documentation/opengrads/doc/udxt/gsudf/index.html">
                                   GsUDF </a> - Writing User Defined Functions as GrADS Scripts
   </li>

   <li><a href="http://opengrads.org/wiki/index.php?title=Gxyat">
                                   GxYAT </a> - A command for creating hight quality images in PNG, Postscript, PDF and SVG formats
   </li>

   <li><a href="Resources/Documentation/opengrads/doc/udxt/libbjt/index.html">
                                   Libbjt </a> - Ben-Jei Tsuang's Function Collection
   </li>

   <li><a href="Resources/Documentation/opengrads/doc/udxt/libmf/index.html">
                                   Libmf </a> - Mike Fiorino's Extension Collection
   </li>

   <li><a href="Resources/Documentation/opengrads/doc/udxt/libipc/index.html">
                                   Libipc </a> - Functions for Inter-process Communication
   </li>
   <li><a href="Resources/Documentation/opengrads/doc/udxt/re/"index.html>
                                   Re </a> - 2D Regridding
   </li>

   <li><a href="Resources/Documentation/opengrads/doc/udxt/saakeskus/index.html">
                                   Saakeskus </a> - Total Totals, Storm-relative Helicity, Potential instability, Wind-chill temperature and other thunderstorm related indices.
   </li>

   <li><a href="Resources/Documentation/opengrads/doc/udxt/shfilt/index.html">
                                   ShFilt </a> - Spherical Harmonic Filtering
   </li>

   </ul>

<h1> <a name="gagui"</a> Using GrADS with Athena Widgets</h1>

GrADS offers a simpler way of creating a Graphical User Interface
(GUI) based on <b>libsx</b>, the <i> Simple X Library</i> by Dominic
Giampaolo. Libsx is a C library layered on top of the Athena widget
set which allows the programming reasonable interfaces with minimum
effort. GrADS includes an interface to <b>libsx</b>, so that users can
enjoy the same simplicity when creating basic graphical user
interfaces in GrADS.

<ul>
<li> <a href=Resources/Documentation/opengrads/doc/gagui/gagui_intro.html>Introduction  </a>
<li> <a href=Resources/Documentation/opengrads/doc/gagui/gagui_run.html>Running the sample script  </a>
<li> <a href=Resources/Documentation/opengrads/doc/gagui/gagui_simple.html>Writing simple scripts  </a>
<li> <a href=Resources/Documentation/opengrads/doc/gagui/gagui_sample.html>Writing the sample script  </a>
<li> <a href=Resources/Documentation/opengrads/doc/gagui/gagui_further.html>Going further...  </a>
<li> <a href=Resources/Documentation/opengrads/doc/gagui/gagui_ref.html>Reference Section  </a>
</ul>


<h1><a name="Bundle" </a> OpenGrADS Bundle Front-ends for Linux/Unix </h1>

The OpenGrADS Bundle provides an unified set of front end scripts for
accessing the basic GrADS application and utilities. The pages below
document the several command line options.

<ul>
   
   <li><a href="Resources/Documentation/opengrads/doc/bundle/grads/index.html">
                grads</a> - the basic grads front end
   </li>

   <li><a href="Resources/Documentation/opengrads/doc/bundle/opengrads/index.html">
                opengrads</a> - grads front end with colorized text and persistent readline history
       
   </li>

   <li><a href="Resources/Documentation/opengrads/doc/bundle/merra/index.html">
                merra </a> - grads front end starting an Athena GUI for visualizing MERRA files; see this <a href="http://cookbooks.opengrads.org/index.php?title=Recipe-016:_Accessing_MERRA_data_with_a_Graphical_User_Interface"> recipe </a> for a description of the<i> MERRA Visualization Tool </i>.
   </li>

</ul>


<h1><a name="Scripts"> </a> Scripts </h1>

   <ul>
   
   <li><a href="Resources/Documentation/opengrads/doc/scripts/lats4d/index.html">
                                  LATS4D</a> - A script for subsetting, reformatting and regriding.
</li>
   </ul>


<? include ("doc/opengrads/footer.php"); ?>
