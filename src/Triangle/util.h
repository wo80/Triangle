/*****************************************************************************/
/*                                                                           */
/*  interpolate()   Interpolate the attibutes of 'newvertex' in the triangle */
/*                  given by 'org', 'dest' and 'apex'.                       */
/*                                                                           */
/*****************************************************************************/

void interpolate(vertex newvertex, vertex org, vertex dest, vertex apex, int nextras)
{
  REAL xdo, ydo, xao, yao;
  REAL denominator;
  REAL dx, dy;
  REAL xi, eta;
  int i;

  xdo = dest[0] - org[0];
  ydo = dest[1] - org[1];
  xao = apex[0] - org[0];
  yao = apex[1] - org[1];

  denominator = 0.5 / (xdo * yao - xao * ydo);

  dx = newvertex[0] - org[0];
  dy = newvertex[1] - org[1];

  // To interpolate vertex attributes for the new vertex, define a
  // coordinate system with a xi-axis directed from the triangle's
  // origin to its destination, and an eta-axis, directed from its
  // origin to its apex.
  xi = (yao * dx - xao * dy) * (2.0 * denominator);
  eta = (xdo * dy - ydo * dx) * (2.0 * denominator);

  for (i = 2; i < 2 + nextras; i++) {
    // Interpolate the vertex attributes.
    newvertex[i] = org[i] + xi * (dest[i] - org[i]) + eta * (apex[i] - org[i]);
  }
}