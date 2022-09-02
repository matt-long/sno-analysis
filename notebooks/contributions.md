# Emails from model output providers 

## MIROC4-ACTM

On Thu, Mar 3, 2022 at 7:46 PM Prabir Patra (via Google Drive) <drive-shares-dm-noreply@google.com> wrote:
prabir@jamstec.go.jp shared a folder
Unknown profile photo
prabir@jamstec.go.jp has invited you to contribute to the following shared folder:

https://drive.google.com/drive/folders/1pJr0WnMhBySiVlhqPROK-hz92a2mRyQi?usp=sharing_eil_m&ts=62217d69

I am uploading MIROC4-ACTM simulations here. Please have a look and let me know if the results look reasonable to you. Thanks and regards, Prabir

MIROC4-ACTM_T42
You have received this email because prabir@jamstec.go.jp shared a file or folder located in Google Drive with you.	Logo for Google Drive

## Jena_TM3

On Mon, May 9, 2022 at 4:00 AM Christian Rödenbeck <Christian.Roedenbeck@bgc-jena.mpg.de> wrote:
Hi Britt,

finally, finally, I managed to do the APO forward runs. See the links to
the results below. I hope that I did everything the same way as last
time. If not, or if anything looks strange (I only looked at a small
sample of output files), please let me know.

For technical reasons (because the runs have formally been done as CO2
runs to be able to use the CO2 obspack input), an initial condition of

    347.0120 ppm

has been added to the output files (but the actual transport runs have
been done from a zero atmosphere).

Cheers
-Christian

http://www.bgc-jena.mpg.de/~christian.roedenbeck/CarboScope/INVERSION/OUTPUT/OUTPUT14.004+tm3_ncep.1986-2020_fg.co2.vGV7.sGV7.FwdFile_BrittStephensAPO/PRI_apojena_output1.vGV7.sGV7.tar.gz

http://www.bgc-jena.mpg.de/~christian.roedenbeck/CarboScope/INVERSION/OUTPUT/OUTPUT14.004+tm3_ncep.1986-2020_fg.co2.vGV7.sGV7.FwdFile_BrittStephensAPO/PRI_co2gridfed_output1.vGV7.sGV7.tar.gz

http://www.bgc-jena.mpg.de/~christian.roedenbeck/CarboScope/INVERSION/OUTPUT/OUTPUT14.004+tm3_ncep.1986-2020_fg.co2.vGV7.sGV7.FwdFile_BrittStephensAPO/PRI_co2oco2mip_output1.vGV7.sGV7.tar.gz

http://www.bgc-jena.mpg.de/~christian.roedenbeck/CarboScope/INVERSION/OUTPUT/OUTPUT14.004+tm3_ncep.1986-2020_fg.co2.vGV7.sGV7.FwdFile_BrittStephensAPO/PRI_co2cesm_output1.vGV7.sGV7.tar.gz

http://www.bgc-jena.mpg.de/~christian.roedenbeck/CarboScope/INVERSION/OUTPUT/OUTPUT14.004+tm3_ncep.1986-2020_fg.co2.vGV7.sGV7.FwdFile_BrittStephensAPO/PRI_co2somffn_output1.vGV7.sGV7.tar.gz

http://www.bgc-jena.mpg.de/~christian.roedenbeck/CarboScope/INVERSION/OUTPUT/OUTPUT14.004+tm3_ncep.1986-2020_fg.co2.vGV7.sGV7.FwdFile_BrittStephensAPO/PRI_n2cesm_output1.vGV7.sGV7.tar.gz

http://www.bgc-jena.mpg.de/~christian.roedenbeck/CarboScope/INVERSION/OUTPUT/OUTPUT14.004+tm3_ncep.1986-2020_fg.co2.vGV7.sGV7.FwdFile_BrittStephensAPO/PRI_n2era_output1.vGV7.sGV7.tar.gz

http://www.bgc-jena.mpg.de/~christian.roedenbeck/CarboScope/INVERSION/OUTPUT/OUTPUT14.004+tm3_ncep.1986-2020_fg.co2.vGV7.sGV7.FwdFile_BrittStephensAPO/PRI_o2gridfed_output1.vGV7.sGV7.tar.gz

http://www.bgc-jena.mpg.de/~christian.roedenbeck/CarboScope/INVERSION/OUTPUT/OUTPUT14.004+tm3_ncep.1986-2020_fg.co2.vGV7.sGV7.FwdFile_BrittStephensAPO/PRI_o2cesm_output1.vGV7.sGV7.tar.gz

http://www.bgc-jena.mpg.de/~christian.roedenbeck/CarboScope/INVERSION/OUTPUT/OUTPUT14.004+tm3_ncep.1986-2020_fg.co2.vGV7.sGV7.FwdFile_BrittStephensAPO/PRI_o2gk01r16_output1.vGV7.sGV7.tar.gz

## CT_TM5

On Mon, Jul 18, 2022 at 4:19 PM Andy Jacobson <andy.jacobson@noaa.gov> wrote:
Guys,

My 137-level APO run 1989-2021 finished. Can I globus you the data? And which data do you want anyway?

-Andy

--
Andy Jacobson
andy.jacobson@noaa.gov

## CAMS_LMDZ

On Thu, Jun 2, 2022 at 8:17 AM Frederic Chevallier <frederic.chevallier@lsce.ipsl.fr> wrote:
Dear Britt, Matt,

Sorry for the slow reaction, but I encountered severe troubles to run my transport model on the supercomputer here further to a software security update. I have prepared the required outputs (none of the encouraged ones) there:

http://dods.lsce.ipsl.fr/invsat/ForML/APO/

Do not hesitate to tell me if you find anything suspicious.

Cheers,

Frederic

## NIES

On Tue, Aug 16, 2022 at 3:44 AM Shamil MAKSYUTOV <shamil@nies.go.jp> wrote:
Hello Britt,

I prepared an draft output to csv file, at obspack/ship locations with NIESTM/Flexpart model. Run starts from zero mixing ratio on jan 1, 2000, units are mol/mol. 

If you can provide a sample in proper csv file form (eg CAMS/Jena), I will rewrite it accordingly

Thanks for your effort

Best

Shamil

On Tue, Aug 23, 2022 at 7:48 PM Shamil MAKSYUTOV <shamil@nies.go.jp> wrote:
Hi Britt,

Thank you for efforts. About producing data for every record in Obspack file. 

I have a technical problem with processing large number of observations, thus I prepared a sparse output for aircraft data (output every n-th observation, where n is selected to keep the total for file under 30000). 

It is not hard to modify my postprocessing program, to output every record, but with a missing value (or nan) where simulation isn't available. Would it be useful?

Note that the vos ship data are in the same file. With site ids: nies, mirai, ham, and observation id is assigned in the order of data in each csv file (row number).

I may also learn adding the data to the Obspack file, seems manageable

Best regards
Shamil 

On Wed, Aug 24, 2022 at 3:40 AM Shamil MAKSYUTOV <shamil@nies.go.jp> wrote:
Hi Britt,

I tried to make netcdf files as in Obspack, see attached. Packing programs for obspack and ships that use my text file are also included

As for additional output, I have 3 hourly output at fixed locations, will send separately

Best regards
Shamil

On Wed, Aug 24, 2022 at 11:08 PM Shamil MAKSYUTOV <shamil@nies.go.jp> wrote:
Hello Britt,

The station output file is available from nies server: fxs.nies.go.jp
user: apofiles
password: HBakNY0r
file name: nies.stations.v1.nc
folder expires: 09/15

Please let me know if anything looks suspicious

Best regards
Shamil