from flask import Flask
from flask import request

from predict_target import predict_target
app = Flask(__name__)


@app.route("/")
def index():
    input_smiles = request.args.get("input_smiles", "")
    if input_smiles:
        prediction = predict_target(input_smiles)
    else:
        prediction = ""
    return (
        """<form action="" method="get">
                <h1>Input SMILES to predict targets</h1>
                <p>Example: Clc1ccc(cc1Cl)C(=O)NCCN1CCC(CC1)n1c2ccccc2[nH]c1=O</p>
                <p> Top 5 predicted targets will be shown</p>
                input_smiles: <input type="text" name="input_smiles">
                <input type="submit" value="Click here for Predicted targets (Uniprot ID)">
            </form>"""
        + "Prediction: "
        + prediction
    )


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=8080, debug=True)













# predict_target("Clc1ccc(cc1Cl)C(=O)NCCN1CCC(CC1)n1c2ccccc2[nH]c1=O")

