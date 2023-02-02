var documenterSearchIndex = {"docs":
[{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"Pages = [\"postprocessing.md\"]","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"","category":"page"},{"location":"postprocessing/#Postprocessing","page":"Postprocessing","title":"Postprocessing","text":"","category":"section"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"When a simulation completes the solver outputs solution data that can be postprocessed to see additional values and/or plots. ","category":"page"},{"location":"postprocessing/#Using-the-simulation-outputs","page":"Postprocessing","title":"Using the simulation outputs","text":"","category":"section"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"When a simulation is executed using t,zm,X,S,Pb,Sb,Lf,sol = BiofilmSolver(p), which is a line of code within each of the Case files solver will produce all the outputs described below","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"t - array of the solution times\nzm - array of the biofilm grid locations at the end of the simulation t[end]\nX=X(t) - array of the particulate concentrations within the tank as a function of time\nS=S(t) - array of the substrate concentrations within the tank as a function of time\nPb=Pb(z) - array of the particulate volume fractions within the biofilm as a function of location within biofilm at the end of the simulation\nSb=Sb(z) - array of the substrate concentrations within the biofilm as a function of location within biofilm at the end of the simulation\nLf=Lf(t) - array of the biofilm thickness as a function of time\nsol - entire solution from the ODE solver - contains the time history of all the variables.  This variable is difficult to parse but is used by the functions described below.","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"The output can be analyzed using Julia commands.  For example, after running Case1.jl (see Run Biofilm.jl), the maximum substrate concentration could be found ","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"julia> maximum(S)\n65.37595010026911","category":"page"},{"location":"postprocessing/#analyzeBiofilm()-Query-Simulation-at-Specified-Times","page":"Postprocessing","title":"analyzeBiofilm() - Query Simulation at Specified Times","text":"","category":"section"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"The analyzeBiofilm(sol,p,t) function provides a simple way to postprocess data.  This function takes in the computed solution sol and parameters p and analyzes the results and the specified time or times t.  For example, after running Case1.jl (see Run Biofilm.jl) you could look at the solution at t=0.25 days using","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"julia> analyzeBiofilm(sol,p,0.25)\nAnalyzing Single Substrate and Particulate Case\n   Time   |      Bug |   Oxygen | min,max(     Bug) | min,max(  Oxygen) |  Lf [μm] \n    0.250 |      102 |     51.7 |     0.08,    0.08 |     29.1,    51.3 |      545","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"which displays the tank and biofilm particulates and substrates and biofilm thickness at the requested time.  Multiple times can be included in the time parameter.  For example to get the values at t=0, 0.25, 0.5, 0.75, 1.0 we can run","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"julia> analyzeBiofilm(sol,p,0:0.25:1)\nAnalyzing Single Substrate and Particulate Case\n   Time   |      Bug |   Oxygen | min,max(     Bug) | min,max(  Oxygen) |  Lf [μm] \n    0.000 |       10 |       10 |     0.08,    0.08 |        0,       0 |       10\n    0.250 |      102 |     51.7 |     0.08,    0.08 |     29.1,    51.3 |      545\n    0.500 |      256 |     2.94 |     0.08,    0.08 |    0.568,    2.87 |      348\n    0.750 |      257 |     2.93 |     0.08,    0.08 |    0.745,    2.87 |      312\n    1.000 |      257 |     2.93 |     0.08,    0.08 |    0.761,    2.87 |      309","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"Adding an optional argument makePlot=true, i.e., analyzeBiofilm(sol,p,0.25,makePlot=true) will produce a plot of the biofilm quantities at the specified time(s)","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"julia> analyzeBiofilm(sol,p,0.25,makePlot=true)\nAnalyzing Single Substrate and Particulate Case\n   Time   |      Bug |   Oxygen | min,max(     Bug) | min,max(  Oxygen) |  Lf [μm] \n    0.250 |      102 |     51.7 |     0.08,    0.08 |     29.1,    51.3 |      545","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"(Image: Plots from analyzeBiofilm)","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"analyzeBiofilm","category":"page"},{"location":"postprocessing/#Biofilm.analyzeBiofilm","page":"Postprocessing","title":"Biofilm.analyzeBiofilm","text":"analyzeBiofilm(sol,p,t)\nanalyzeBiofilm(sol,p,t,makePlot=true)\n\nTake solution from biofilm solver and outputs variabes and a plot of biofilm variables.\n\n\n\n\n\n","category":"function"},{"location":"postprocessing/#movieBiofilm()-Make-movie-of-biofilm-quantities","page":"Postprocessing","title":"movieBiofilm() - Make movie of biofilm quantities","text":"","category":"section"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"The particulate and substrate concentrations change throughout the simulation, and it is often useful to make movies of how these change over time.  ","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"The movieBiofilm(times) function provides a convenient way to make these movies. For example, to post process Case 5 - Phototroph, which has a light that turns on and off throughout each day we could make a movie of the biofilm conditions during the day and during night.  To make a movie of every 5 days when the light is on, that is at times=1,5,10,... we can use","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"julia> movieBiofilm(sol,p,1:5:t[end],filename=\"phototroph_day.gif\",fps=5)","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"which will analyze each day and combine the results into the movie: (Image: phototroph during day)","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"To produce the results every 5th night, that is at times=0.5, 5.5, ... we can use","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"julia> movieBiofilm(sol,p,0.5:5:t[end],filename=\"phototroph_night.gif\",fps=5)","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"(Image: phototroph during day)","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"Note that during the day the growth rate is zero and the oxygen concentration is 8.6, which is the inflow concentration Sin.  During the day however, the phototroph growth rate is non-zero, and the growth produces oxygen. ","category":"page"},{"location":"postprocessing/","page":"Postprocessing","title":"Postprocessing","text":"movieBiofilm","category":"page"},{"location":"postprocessing/#Biofilm.movieBiofilm","page":"Postprocessing","title":"Biofilm.movieBiofilm","text":"movieBiofilm(sol,p,times)\nmovieBiofilm(sol,p,times,filename=\"anim.gif\", fps=20)\n\nMake a movie of the biofilm particulate volume fraction, substrate concentration, and particulate growthrates at the specified times.\n\nOptional arguments \n\nfilename: name and type of output, i.e., \"biofilm.mp4\", \"biofilm.gif\"\nframerate\n\nExamples:\n\nCreate movie with t=0,1,...,10\n\njulia> movieBiofilm(sol,p,0:1:10)\n\nCreate movie with specified filename and framerate\n\njulia> movieBiofilm(sol,p,0:1:10,filename=\"biofilm.gif\",fps=10)\n\n\n\n\n\n","category":"function"},{"location":"theory/#Theory","page":"Solver Theory","title":"Theory","text":"","category":"section"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"Biofilm.jl simulates a one-dimensional biofilm within a stirred tank reactor.  The dependent variables include the ","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"tank particulates (biomass) concentration(s) X,\ntank substrate concentrations S,\nbiofilm particulate (biomass) volume fractions P_b,\nbiofilm substrate concentrations S_b, and\nbiofilm thickness L_f. ","category":"page"},{"location":"theory/#Tank-Equations","page":"Solver Theory","title":"Tank Equations","text":"","category":"section"},{"location":"theory/#Particulates","page":"Solver Theory","title":"Particulates","text":"","category":"section"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"The governing equation describing the particulate concentrations in the tank environment is","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"fracd X_tjdt = mu_j(mathbfS_t) X_tj - fracQ X_tjV + fracv_mathrmdet A X_bj(L_f)V + mathrmsrc_Xj","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"for j=1dotsN_x, where t is time, mu_j(mathbfS_t) is the growthrate of the j^mathrmth particulate, Q is the flowrate, V is the volume of the tank, v_mathrmdet=K_mathrmdet L_f^2 is the detachment velocity, A is the area of the biofilm, X_bj(L_f) is the j^mathrmth particulate concentration at the top of the biofilm, and mathrmsrc_Xj is the source term for the j^mathrmth particulate. ","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"The terms on the right-hand-side (RHS) are ","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"the growth of the particulate in the tank, \ntransport due to flow out of the tank, \ntransfer of particulates from the biofilm to the tank due to detachment, and\nsource term.","category":"page"},{"location":"theory/#Substrates","page":"Solver Theory","title":"Substrates","text":"","category":"section"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"The governing equation describing the substrate concentrations in the tank environment is","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"fracd S_tkdt = -sum_j=1^N_x fracmu_j(mathbfS_t) X_tjY_jk + fracQ S_mathrminkV - fracQ S_tkV + fracA S_mathrmfluxkV + mathrmsrc_Sk","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"for k=1dotsN_s, where S_mathrmfluxk is the flux of substrates from the biofilm into the tank, and mathrmsrc_Sk is the source term for the k^mathrmth substrate. ","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"The terms on the right-hand-side (RHS) are ","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"consumption of substrates due to the growth of the particulate in the tank, \ntransport due to flow into the tank, \ntransport due to flow out of the tank,\ntransfer of substrates into the biofilm due to diffusion, and\nsource term.","category":"page"},{"location":"theory/#Biofilm-Equations","page":"Solver Theory","title":"Biofilm Equations","text":"","category":"section"},{"location":"theory/#Particulates-2","page":"Solver Theory","title":"Particulates","text":"","category":"section"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"The governing equations describing the biofilm environment are","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"fracd P_bjidt = \nmu_j(mathbfS_bi) P_bji \n- fracd v_i P_bjidz \n+ fracmathrmsrc_Xjirho_j","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"for j=1dotsN_x and i=1dotsN_z. Where P_bji is the j^mathrmth particulate at the i^mathrmth grid point within the biofilm. ","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"The terms on the right-hand-side (RHS) are ","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"the growth of the particulate in the biofilm, \ntransport through the biofilm due to the growth velocity v_i, and \nsource term of particulate at i^mathrmth location in biofilm.","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"The growth velocity v_i is the rate of flow through the biofilm due to growth deeper within the biofilm and is defined with","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"v_i=  int_z=0^z_isum_j=1^N_x frac1P_mathrmtotleft(mu_j(mathbfS_bi) P_bji + fracmathrmsrc_Xjrho_jright) dz","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"where P_mathrmtot=sum_j=1^N_xP_bj","category":"page"},{"location":"theory/#Substrates-2","page":"Solver Theory","title":"Substrates","text":"","category":"section"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"fracd S_bkidt = \nD_ekfracd^2 S_bkidz^2 \n- sum_j=1^N_x fracmu_j(mathbfS_bi) X_bjiY_jk\n+ mathrmsrc_Ski","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"for k=1dotsN_s and i=1dotsN_z.","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"The terms on the right-hand-side (RHS) are ","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"diffusion of substrates in the biofilm,\nconsumption of substrates due to the growth of the particulate in the biofilm, and\nsource term of substrate at i^mathrmth location in biofilm.","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"The diffusion term with a second derivative w.r.t. z requires boundary conditions at the top and bottom of the biofilm.  A zero-flux (zero first-derivative) condition is used at the bottom of the biofilm.  At the top of the biofilm the diffusion through the boundary layer is matched with the diffusion into the biofilm, i.e.,","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"D_aqkfracd^2 S_kdz^2 = D_ekfracd^2 S_bkN_zdz^2 ","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"for k=1dotsN_s, where D_e is the diffusion coefficient in the biofilm and D_mathrmaq is the diffusion coefficient in the boundary layer. ","category":"page"},{"location":"theory/#Biofilm-Thickness","page":"Solver Theory","title":"Biofilm Thickness","text":"","category":"section"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"The thickness of the biofilm L_f is described by ","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"fracd L_fdt = v_N_z - v_mathrmdet","category":"page"},{"location":"theory/","page":"Solver Theory","title":"Solver Theory","text":"where the first term on the RHS is the growth velocity at the top of the biofilm (see Biofilm Particulates) and the second term is the detachment velocity modeled with v_mathrmdet=K_mathrmdet L_f^2","category":"page"},{"location":"examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Pages = [\"examples.md\"]","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"","category":"page"},{"location":"examples/#Case-1-Single-Substrate-and-Particulate-Case","page":"Examples","title":"Case 1 - Single Substrate and Particulate Case","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Download link: Case1.jl (right-click & Save Link As)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"This case is the same the case shown on the Installation page.  It simulates a single particulate (\"Bug\") and a single substrate (\"Oxygen\").  The particulate has a growthrate of ","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"mu = mu_mathrmmax fracSK_M + S","category":"page"},{"location":"examples/#Case-2-Multiple-Substrates","page":"Examples","title":"Case 2 - Multiple Substrates","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Download link: Case2.jl (right-click & Save Link As)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"This is a simple example of how multiple substrates can be simulated.  Substrate 1 is used by \"Bug\" while substrate 2 is completely indpendent (not very interesting).","category":"page"},{"location":"examples/#Case-3-Live/Dead-Bugs","page":"Examples","title":"Case 3 - Live/Dead Bugs","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Download link: Case3.jl (right-click & Save Link As)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"This example has living bugs that use the substrate to grow.  The bugs die and the concentration of dead bugs is also computed. The source term, src, is used to transfer living bugs to dead bugs.","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"mathrmsrc_X_1 = mathrmsrc_mathrmliving = -b X_1\nmathrmsrc_X_2 = mathrmsrc_mathrmdead = b X_1","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"where b=0.1 is a constant.  The first source term reduces living bugs and the second adds these to the dead bugs. ","category":"page"},{"location":"examples/#Case-4-Multiple-Particulate-and-Substrates","page":"Examples","title":"Case 4 - Multiple Particulate and Substrates","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Download link: Case4.jl (right-click & Save Link As)","category":"page"},{"location":"examples/#Case-5-Phototroph","page":"Examples","title":"Case 5 - Phototroph","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Download link: Case5_Photroph.jl (right-click & Save Link As)","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"This case shows how you can simulate a species that grows in light.  The light turns off and on throughout each day and also only penetrates the top of the biofilm.  Since the light has discontinuities, the case has an additional parameter","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"    discontinuityPeriod=0.25,","category":"page"},{"location":"examples/","page":"Examples","title":"Examples","text":"that tells the solver the discontinuities can occur every 0.25 days. ","category":"page"},{"location":"examples/#Case-6-Three-Species","page":"Examples","title":"Case 6 - Three Species","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Download link: Case6_threeSpecies.jl (right-click & Save Link As)","category":"page"},{"location":"examples/#Case-7-SOB-and-SRB","page":"Examples","title":"Case 7 - SOB & SRB","text":"","category":"section"},{"location":"examples/","page":"Examples","title":"Examples","text":"Download link: Case7_SOB_SRB.jl (right-click & Save Link As)","category":"page"},{"location":"installation/#Installation","page":"Installation","title":"Installation","text":"","category":"section"},{"location":"installation/#Download-and-Install-Julia","page":"Installation","title":"Download & Install Julia","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"Download the latest version of Julia from julialang.org\nInstall Julia following the provided directions on the help page.  The default installation directory should be fine.  On Windows be sure to select the option to add Julia to the path to allow VS Code (next section) to find Julia.  On Mac or Linux be sure to follow the instructions on adding Julia to PATH.\nLaunch Julia by finding the icon where you installed the program.  If you have been successful, you should see the REPL like this image (version number may be different)","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"(Image: Julia REPL)","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"tip: Tip\nTry running a simple code like julia> 5+5 to see that you can use Julia.  You could even try learning more about the language by doing a tutorial like this one: \"From zero to Julia!\". ","category":"page"},{"location":"installation/#VS-Code-GUI-for-Julia","page":"Installation","title":"VS Code - GUI for Julia","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"Julia is a programming languages and you can use it in many ways.  One popular way to run Julia code is through Visual Studio Code, which allows you to edit, run, and see results from a code.","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"Download, install, and open VS Code.  You can find the installers on code.visualstudio.com\nInstall the Julia extension in VS Code by selecting File [Code on Mac] >> Preferences >> Extensions, then search for Julia.  You should see an extension called Julia provided by julialang.  Click the Install button. (Image: Julia Extension)\nStart the Julia REPL in VS Code by opening the command-pallet using Ctrl-Shift-p on Windows or Cmd-Shfit-p on Mac then searching for and running Julia: Start REPL.  This will start Julia inside of VS Code.","category":"page"},{"location":"installation/#Add-the-Biofilm.jl-package","page":"Installation","title":"Add the Biofilm.jl package","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"Run the following command in the Julia REPL within VS Code.\njulia> using Pkg; Pkg.add(url=\"https://github.com/markowkes/Biofilm.jl\")","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"This will download Biofilm.jl and all the dependencies and can take several minutes to finish.","category":"page"},{"location":"installation/#Run-Biofilm.jl","page":"Installation","title":"Run Biofilm.jl","text":"","category":"section"},{"location":"installation/","page":"Installation","title":"Installation","text":"Try running the examples.","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"Download Case1.jl by right clicking and choosing Save Link As. Save the file to your Downloads folder (or other location of your choice)\nRun the case by opening the file in VS Code and clicking the play button in the top right corner.  You should see output to the REPL and a plot like this: (Image: Case 1 Run Button)","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"The top row shows the biomass (particulate) and substrate concentrations and the biofilm thickness as a function of time.  The bottom row shows the particulate volume fraction and substrate concentrations as a function of position within the biofilm at the output time (t=1.0 as indicated by the title).","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"Try editing the case file or explore other Examples","category":"page"},{"location":"installation/","page":"Installation","title":"Installation","text":"note: Note\nThe first time you run Biofilm.jl it will require time to compile. Subsequent runs should occur much faster. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"This package models the dynamics of a biofilm using the Julia programming language.  ","category":"page"},{"location":"","page":"Home","title":"Home","text":" Pages = [\n   \"index.md\",\n   \"installation.md\",\n   \"parameters.md\",\n   \"examples.md\",\n   \"postprocessing.md\",\n   \"theory.md\",\n ]","category":"page"},{"location":"parameters/#Case-Parameters","page":"Case Parameters","title":"Case Parameters","text":"","category":"section"},{"location":"parameters/","page":"Case Parameters","title":"Case Parameters","text":"Different biofilms can be modeled by modify the Case1.jl file.  This is a text file and can be opened and edited using Notepad (Windows), textEdit (Mac), or similar. ","category":"page"},{"location":"parameters/","page":"Case Parameters","title":"Case Parameters","text":"The Case file defines a param struct with the following fields","category":"page"},{"location":"parameters/#Simulation-Parameters","page":"Case Parameters","title":"Simulation Parameters","text":"","category":"section"},{"location":"parameters/","page":"Case Parameters","title":"Case Parameters","text":"Title - Description of the case, used, e.g., on the title of plots \ntFinal - Simulation is performed from t=0 to t=t_mathrmFinal [days]\ntol - Tolerance used for differential equation solver.  Solution will have error le tol.\noutPeriod - Period in days between outputs during simulation.  Note that the solver will take smaller timesteps then outPeriod and the entire solution will be available when the solver completes (that is outPeriod has no impact on the solution). [days]\nplotPeriod - [default=outPeriod] - Period in days between plot renderings during the simulation.  For long simulations you can increase plotPeriod to reduce the number of plots that are created and speedup the simulation. Note: plotPeriod is required to be a multiple of outPeriod. [days]","category":"page"},{"location":"parameters/#Particulate-Parameters","page":"Case Parameters","title":"Particulate Parameters","text":"","category":"section"},{"location":"parameters/","page":"Case Parameters","title":"Case Parameters","text":"XNames - Names of the particulates used in plots and other outputs\nXto - Initial condition of particulate concentration(s) in tank environment [g/m^3]\nPbo - Initial condition of particulate volume fraction(s) in biofilm [-]\nrho - Particulate densities [g/m^3]\nKdet - Particulate detachment coefficient [1/(m days)]\nsrc - Array of source terms for each particulate.  Each source term can be a function of substrate or particulate concentrations.  See Case 4 for an example of using src. [g/(m^3 s)]\nmu - Growthrates for each particulate.  Each growthrate can be a function of substrate and particulate concentrations, biofilm thickness, time, position within biofilm, and parameters.  See example for additional details. [1/days]\ndiscontinuityPeriod - [Optional] - Period in days between discontinuities in mu or src.  If mu or src are continuous do not set this optional parameter. [days]","category":"page"},{"location":"parameters/#Substrate-Parameters","page":"Case Parameters","title":"Substrate Parameters","text":"","category":"section"},{"location":"parameters/","page":"Case Parameters","title":"Case Parameters","text":"SNames - Names of the substrates used in plots and other outputs\nSo - Initial condition of substrate concentration(s) in tank environment [g/m^3]\nSbo - Initial condition of substrate concentrations(s) in biofilm [g/m^3]\nYxs - Array of biomass yield coefficients [g_mathrmX/g_mathrmS]\nY_xs = beginbmatrix\nfracDelta X_1Delta S_1  fracDelta X_1Delta S_2  cdots  fracDelta X_1Delta S_N05em\nfracDelta X_2Delta S_1  fracDelta X_2Delta S_2  cdots  fracDelta X_2Delta S_N05em\nvdots  vdots  ddots  vdots 05em\nfracDelta X_MDelta S_1  fracDelta X_MDelta S_2  cdots  fracDelta X_MDelta S_N\nendbmatrix\nfor N substrates and M particulates.  If X_i does not depend on S_j, set displaystyle fracDelta X_iDelta S_j=0.\nDaq - Diffusion coefficient of substrates through a boundary at the top of the biofilm.  If you do not want a boundary layer set LL=0, then these values will not impact the solution. [m^2/day]\nDe - Diffusion coefficient of substrates through a biofilm [m^2/day]","category":"page"},{"location":"parameters/#Tank-Parameters","page":"Case Parameters","title":"Tank Parameters","text":"","category":"section"},{"location":"parameters/","page":"Case Parameters","title":"Case Parameters","text":"V - Volume of tank [m^3]\nS - Surface area of biofilm [m^2]\nQ - Flowrate through tank [m^3/day]","category":"page"},{"location":"parameters/#Biofilm-Parameters","page":"Case Parameters","title":"Biofilm Parameters","text":"","category":"section"},{"location":"parameters/","page":"Case Parameters","title":"Case Parameters","text":"Nz - Number of grid points used to discretize the biofilm\nLfo - Initial condition for biofilm thickness [m]\nLL - Thickness of boundary layer at top of biofilm [m]","category":"page"}]
}
