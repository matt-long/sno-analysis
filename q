diff --git a/notebooks/_cmip6-sno-precompute.ipynb b/notebooks/_cmip6-sno-precompute.ipynb
index 7a5db5e..817373d 100644
--- a/notebooks/_cmip6-sno-precompute.ipynb
+++ b/notebooks/_cmip6-sno-precompute.ipynb
@@ -58,7 +58,7 @@
     {
      "data": {
       "text/html": [
-       "<p><strong>glade-cmip6 catalog with 5226 dataset(s) from 2303258 asset(s)</strong>:</p> <div>\n",
+       "<p><strong>glade-cmip6 catalog with 5260 dataset(s) from 2325959 asset(s)</strong>:</p> <div>\n",
        "<style scoped>\n",
        "    .dataframe tbody tr th:only-of-type {\n",
        "        vertical-align: middle;\n",
@@ -118,15 +118,15 @@
        "    </tr>\n",
        "    <tr>\n",
        "      <th>version</th>\n",
-       "      <td>595</td>\n",
+       "      <td>597</td>\n",
        "    </tr>\n",
        "    <tr>\n",
        "      <th>time_range</th>\n",
-       "      <td>39054</td>\n",
+       "      <td>44033</td>\n",
        "    </tr>\n",
        "    <tr>\n",
        "      <th>path</th>\n",
-       "      <td>2303258</td>\n",
+       "      <td>2325959</td>\n",
        "    </tr>\n",
        "  </tbody>\n",
        "</table>\n",
@@ -154,8 +154,8 @@
      "name": "stdout",
      "output_type": "stream",
      "text": [
-      "CPU times: user 6.24 s, sys: 340 ms, total: 6.58 s\n",
-      "Wall time: 6.71 s\n"
+      "CPU times: user 6.86 s, sys: 536 ms, total: 7.4 s\n",
+      "Wall time: 7.39 s\n"
      ]
     },
     {
@@ -285,7 +285,7 @@
        "      <td>...</td>\n",
        "    </tr>\n",
        "    <tr>\n",
-       "      <th>2303253</th>\n",
+       "      <th>2325954</th>\n",
        "      <td>ScenarioMIP</td>\n",
        "      <td>MRI</td>\n",
        "      <td>MRI-ESM2-0</td>\n",
@@ -300,7 +300,7 @@
        "      <td>/glade/collections/cmip/CMIP6/gcm/ScenarioMIP/...</td>\n",
        "    </tr>\n",
        "    <tr>\n",
-       "      <th>2303254</th>\n",
+       "      <th>2325955</th>\n",
        "      <td>ScenarioMIP</td>\n",
        "      <td>MRI</td>\n",
        "      <td>MRI-ESM2-0</td>\n",
@@ -315,7 +315,7 @@
        "      <td>/glade/collections/cmip/CMIP6/gcm/ScenarioMIP/...</td>\n",
        "    </tr>\n",
        "    <tr>\n",
-       "      <th>2303255</th>\n",
+       "      <th>2325956</th>\n",
        "      <td>ScenarioMIP</td>\n",
        "      <td>MRI</td>\n",
        "      <td>MRI-ESM2-0</td>\n",
@@ -330,7 +330,7 @@
        "      <td>/glade/collections/cmip/CMIP6/gcm/ScenarioMIP/...</td>\n",
        "    </tr>\n",
        "    <tr>\n",
-       "      <th>2303256</th>\n",
+       "      <th>2325957</th>\n",
        "      <td>ScenarioMIP</td>\n",
        "      <td>MRI</td>\n",
        "      <td>MRI-ESM2-0</td>\n",
@@ -345,7 +345,7 @@
        "      <td>/glade/collections/cmip/CMIP6/gcm/ScenarioMIP/...</td>\n",
        "    </tr>\n",
        "    <tr>\n",
-       "      <th>2303257</th>\n",
+       "      <th>2325958</th>\n",
        "      <td>ScenarioMIP</td>\n",
        "      <td>MRI</td>\n",
        "      <td>MRI-ESM2-0</td>\n",
@@ -361,7 +361,7 @@
        "    </tr>\n",
        "  </tbody>\n",
        "</table>\n",
-       "<p>2303258 rows × 12 columns</p>\n",
+       "<p>2325959 rows × 12 columns</p>\n",
        "</div>"
       ],
       "text/plain": [
@@ -372,11 +372,11 @@
        "3         AerChemMIP            BCC    BCC-ESM1  ssp370-lowNTCF  r2i1p1f1   \n",
        "4         AerChemMIP            BCC    BCC-ESM1  ssp370-lowNTCF  r3i1p1f1   \n",
        "...              ...            ...         ...             ...       ...   \n",
-       "2303253  ScenarioMIP            MRI  MRI-ESM2-0          ssp585  r1i1p1f1   \n",
-       "2303254  ScenarioMIP            MRI  MRI-ESM2-0          ssp585  r1i1p1f1   \n",
-       "2303255  ScenarioMIP            MRI  MRI-ESM2-0          ssp585  r1i1p1f1   \n",
-       "2303256  ScenarioMIP            MRI  MRI-ESM2-0          ssp585  r1i1p1f1   \n",
-       "2303257  ScenarioMIP            MRI  MRI-ESM2-0          ssp585  r1i1p1f1   \n",
+       "2325954  ScenarioMIP            MRI  MRI-ESM2-0          ssp585  r1i1p1f1   \n",
+       "2325955  ScenarioMIP            MRI  MRI-ESM2-0          ssp585  r1i1p1f1   \n",
+       "2325956  ScenarioMIP            MRI  MRI-ESM2-0          ssp585  r1i1p1f1   \n",
+       "2325957  ScenarioMIP            MRI  MRI-ESM2-0          ssp585  r1i1p1f1   \n",
+       "2325958  ScenarioMIP            MRI  MRI-ESM2-0          ssp585  r1i1p1f1   \n",
        "\n",
        "        table_id variable_id grid_label  dcpp_init_year    version  \\\n",
        "0            day        rsds         gn             NaN  v20190624   \n",
@@ -385,11 +385,11 @@
        "3            day      tasmax         gn             NaN  v20190624   \n",
        "4            day        rsds         gn             NaN  v20190624   \n",
        "...          ...         ...        ...             ...        ...   \n",
-       "2303253   6hrLev          va         gn             NaN  v20190625   \n",
-       "2303254   6hrLev          va         gn             NaN  v20190625   \n",
-       "2303255     Oday         tos         gn             NaN  v20210329   \n",
-       "2303256     Oday         tos         gn             NaN  v20210329   \n",
-       "2303257       fx        orog         gn             NaN  v20190603   \n",
+       "2325954   6hrLev          va         gn             NaN  v20190625   \n",
+       "2325955   6hrLev          va         gn             NaN  v20190625   \n",
+       "2325956     Oday         tos         gn             NaN  v20210329   \n",
+       "2325957     Oday         tos         gn             NaN  v20210329   \n",
+       "2325958       fx        orog         gn             NaN  v20190603   \n",
        "\n",
        "                        time_range  \\\n",
        "0                20150101-20551231   \n",
@@ -398,11 +398,11 @@
        "3                20150101-20551231   \n",
        "4                20150101-20551231   \n",
        "...                            ...   \n",
-       "2303253  209803010000-209803311800   \n",
-       "2303254  209804010000-209804301800   \n",
-       "2303255          20150101-20641231   \n",
-       "2303256          20650101-21001231   \n",
-       "2303257                        NaN   \n",
+       "2325954  209803010000-209803311800   \n",
+       "2325955  209804010000-209804301800   \n",
+       "2325956          20150101-20641231   \n",
+       "2325957          20650101-21001231   \n",
+       "2325958                        NaN   \n",
        "\n",
        "                                                      path  \n",
        "0        /glade/collections/cmip/CMIP6/AerChemMIP/BCC/B...  \n",
@@ -411,13 +411,13 @@
        "3        /glade/collections/cmip/CMIP6/AerChemMIP/BCC/B...  \n",
        "4        /glade/collections/cmip/CMIP6/AerChemMIP/BCC/B...  \n",
        "...                                                    ...  \n",
-       "2303253  /glade/collections/cmip/CMIP6/gcm/ScenarioMIP/...  \n",
-       "2303254  /glade/collections/cmip/CMIP6/gcm/ScenarioMIP/...  \n",
-       "2303255  /glade/collections/cmip/CMIP6/gcm/ScenarioMIP/...  \n",
-       "2303256  /glade/collections/cmip/CMIP6/gcm/ScenarioMIP/...  \n",
-       "2303257  /glade/collections/cmip/CMIP6/gcm/ScenarioMIP/...  \n",
+       "2325954  /glade/collections/cmip/CMIP6/gcm/ScenarioMIP/...  \n",
+       "2325955  /glade/collections/cmip/CMIP6/gcm/ScenarioMIP/...  \n",
+       "2325956  /glade/collections/cmip/CMIP6/gcm/ScenarioMIP/...  \n",
+       "2325957  /glade/collections/cmip/CMIP6/gcm/ScenarioMIP/...  \n",
+       "2325958  /glade/collections/cmip/CMIP6/gcm/ScenarioMIP/...  \n",
        "\n",
-       "[2303258 rows x 12 columns]"
+       "[2325959 rows x 12 columns]"
       ]
      },
      "execution_count": 4,
@@ -548,7 +548,7 @@
      "name": "stderr",
      "output_type": "stream",
      "text": [
-      "/glade/work/mclong/miniconda3/envs/sno/lib/python3.7/site-packages/xarray/conventions.py:520: SerializationWarning: variable 'areacello' has multiple fill values {1e+20, 1e+20}, decoding all values to NaN.\n",
+      "/glade/u/home/stephens/miniconda3/envs/sno/lib/python3.7/site-packages/xarray/conventions.py:520: SerializationWarning: variable 'areacello' has multiple fill values {1e+20, 1e+20}, decoding all values to NaN.\n",
       "  decode_timedelta=decode_timedelta,\n"
      ]
     },
@@ -3293,7 +3293,7 @@
    "name": "python",
    "nbconvert_exporter": "python",
    "pygments_lexer": "ipython3",
-   "version": "3.7.11"
+   "version": "3.7.10"
   }
  },
  "nbformat": 4,
