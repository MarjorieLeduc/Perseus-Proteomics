
using System;
using System.Collections.Generic;
using BaseLibS.Graph;
using BaseLibS.Param;
using PerseusApi.Document;
using PerseusApi.Generic;
using PerseusApi.Matrix;

namespace Plugin3P5_CombineMultipleCategoricalColumns
{
	public class HeadProcessing : IMatrixProcessing
	{
		public bool HasButton => false;
		public Bitmap2 DisplayImage => null;
		public string Description => "Combine the terms in two categorical columns to form a new categorical column.";
		public string HelpOutput => "A new categorical column is generated with combined terms.";
		public string[] HelpSupplTables => new string[0];
		public int NumSupplTables => 0;
		public string Name => "Combine multiple categorical columns 3P5";
		public string Heading => "Rearrange";
		public bool IsActive => true;
		public float DisplayRank => 17.5f;
		public string[] HelpDocuments => new string[0];
		public int NumDocuments => 0;

		public string Url
			=> "http://coxdocs.org/doku.php?id=perseus:user:activities:MatrixProcessing:Rearrange:CombineCategoricalColumns";

		public int GetMaxThreads(Parameters parameters)
		{
			return 1;
		}

		public Parameters GetParameters(IMatrixData mdata, ref string errorString)
		{
			return
				new Parameters(new SingleChoiceParam("First column", 0) { Values = mdata.CategoryColumnNames },
				new SingleChoiceParam("Second column", 0) { Values = mdata.CategoryColumnNames },
				new BoolWithSubParams("Third column", false)
				{
					Help = "Specify if there is a third column to concat.",
					SubParamsFalse = new Parameters(new Parameter[] { }),
					SubParamsTrue =
						new Parameters(new Parameter[]
							{new SingleChoiceParam("Name of the third column",0){Values = mdata.CategoryColumnNames}})
				},

				new BoolWithSubParams("Fourth column", false)
				{
					Help = "Specify if there is a fourth column to concat.",
					SubParamsFalse = new Parameters(new Parameter[] { }),
					SubParamsTrue =
						new Parameters(new Parameter[]
							{new SingleChoiceParam("Name of the fourth column",0){Values = mdata.CategoryColumnNames}})
				},

				new BoolWithSubParams("Fifth column", false)
				{
					Help = "Specify if there is a fifth column to concat.",
					SubParamsFalse = new Parameters(new Parameter[] { }),
					SubParamsTrue =
						new Parameters(new Parameter[]
							{new SingleChoiceParam("Name of the fifth column",0){Values = mdata.CategoryColumnNames}})
				});
		}

		public void ProcessData(IMatrixData mdata, Parameters param, ref IMatrixData[] supplTables,
			ref IDocumentData[] documents, ProcessInfo processInfo)
		{
			if (mdata.CategoryColumnCount < 2)
			{
				processInfo.ErrString = "There are less than two categorical columns available.";
				return;
			}

			// combine colonne 1 + 2
			int colInd1 = param.GetParam<int>("First column").Value;
			int colInd2 = param.GetParam<int>("Second column").Value;
			string[][] col1 = mdata.GetCategoryColumnAt(colInd1);
			string[][] col2 = mdata.GetCategoryColumnAt(colInd2);
			string[][] result = new string[col1.Length][];
			for (int i = 0; i < result.Length; i++)
			{
				result[i] = CombineTerms(col1[i], col2[i]);
			}
			string colName = mdata.CategoryColumnNames[colInd1] + "_" + mdata.CategoryColumnNames[colInd2];


			// combine colonne1_2 + 3
			if (param.GetParamWithSubParams<bool>("Third column").Value == true)
			{
				colInd2 = param.GetParamWithSubParams<bool>("Third column").GetSubParameters().GetParam<int>("Name of the third column").Value;
				col1 = result;
				col2 = mdata.GetCategoryColumnAt(colInd2);
				for (int i = 0; i < result.Length; i++)
				{
					result[i] = CombineTerms(col1[i], col2[i]);
				}
				colName = colName + "_" + mdata.CategoryColumnNames[colInd2];
			}


			// combine colonne1_2_3 + 4
			if (param.GetParamWithSubParams<bool>("Fourth column").Value == true)
			{
				colInd2 = param.GetParamWithSubParams<bool>("Fourth column").GetSubParameters().GetParam<int>("Name of the fourth column").Value;
				col1 = result;
				col2 = mdata.GetCategoryColumnAt(colInd2);
				for (int i = 0; i < result.Length; i++)
				{
					result[i] = CombineTerms(col1[i], col2[i]);
				}
				colName = colName + "_" + mdata.CategoryColumnNames[colInd2];
			}


			// combine colonne1_2_3_4 + 5
			if (param.GetParamWithSubParams<bool>("Fifth column").Value == true)
			{
				colInd2 = param.GetParamWithSubParams<bool>("Fifth column").GetSubParameters().GetParam<int>("Name of the fifth column").Value;
				col1 = result;
				col2 = mdata.GetCategoryColumnAt(colInd2);
				for (int i = 0; i < result.Length; i++)
				{
					result[i] = CombineTerms(col1[i], col2[i]);
				}
				colName = colName + "_" + mdata.CategoryColumnNames[colInd2];
			}


			mdata.AddCategoryColumn(colName, "", result);

		}

		private static string[] CombineTerms(ICollection<string> x, ICollection<string> y)
		{
			string[] result = new string[x.Count * y.Count];
			int count = 0;
			foreach (string t in x)
			{
				foreach (string t1 in y)
				{
					result[count++] = t + "_" + t1;
				}
			}
			Array.Sort(result);
			return result;
		}
	}
}


