/*      POA_ABS     absolurte height and orientation
 *      POA_ZREL    height is relative to surface level
 *      POA_SURFREL height and orientation relative to surface level
 */
typedef enum POA_SETTINGS {POA_ABS, POA_ZREL, POA_SURFREL} POA_SETTINGS; 


typedef struct simulation_config {
	AOI_Model_Data M; // angle of incidence effects
	POA_SETTINGS P; // setup plane of array
} simulation_config;

simulation_config InitConf();
void InitVars();
void ClearVars();
int lookupvar(char *name);
int LookupTopology(char *name, topology **T);
int LookupSkyDome(char *name, sky_grid **S);
int LookupVec(char *name, sky_pos **v);
int LookupSimConf(char *name, simulation_config **C);
int AddTopology(char *name, topology T);
int AddSkyDome(char *name, sky_grid S);
int AddVec(char *name, sky_pos v);
int AddSimConf(char *name, simulation_config C);
int RMVar(char *name);

void ListVars();
