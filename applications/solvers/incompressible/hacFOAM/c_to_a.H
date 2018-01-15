/*---------------------------------------------------------------------------*\
c_to_a.H

compute average velocity for each cell in CFD and send values to MD
\*---------------------------------------------------------------------------*/

scalarField scounter(nC2ADataPts_p1*nInterfaces, 0);
scalarField vx_sum(nC2ADataPts_p1*nInterfaces, 0);
scalarField vy_sum(nC2ADataPts_p1*nInterfaces, 0);
scalarField vz_sum(nC2ADataPts_p1*nInterfaces, 0);

for (int i=0; i<nInterfaces; i++)
{
  for (int j=0; j<nC2ADataPts_p1; j++)
  {
    for(int k=0; k< C2A_Cells[j].size(); k++)
    {
      int ind=C2A_Cells[j][k];

      int ij=i*nC2ADataPts_p1+j;
      vx_sum[ij]=vx_sum[ij]+U[ind].component(0);
      vy_sum[ij]=vy_sum[ij]+U[ind].component(1);
      vz_sum[ij]=vz_sum[ij]+U[ind].component(2);
      scounter[ij]=scounter[ij]+1;
    }
  }
}

reduce(scounter, sumOp<scalarField>());
reduce(vx_sum, sumOp<scalarField>());
reduce(vy_sum, sumOp<scalarField>());
reduce(vz_sum, sumOp<scalarField>());

scalarField vx_per_bin = vx_sum/scounter;
scalarField vy_per_bin = vy_sum/scounter;
scalarField vz_per_bin = vz_sum/scounter;

for (int i=0; i<nInterfaces; i++) 
{		
  int ifsIndex=i+1;
  c_2_a_velMap_x[ifsIndex]=vx_per_bin;
  c_2_a_velMap_y[ifsIndex]=vy_per_bin;
}
