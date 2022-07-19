# -*- coding: utf-8 -*-
"""hw3_2.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1oI-nnXz99BrPNHZzMKYlveRaZRr_xMN5
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim

from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA

"""# Autoencoder
In this question, we will try to implement **autoencoder(AE)** for dimensionality reduction using `pytorch`. We will test our model on peripheral blood mononuclear cells(PBMCs) single-cell gene expression dataset.

An autoencoder consists of two parts: 

* The encoder network encodes the high dimensional data into low dimensional embedding. 
* The decoder network reconstructs the input high dimensional data from the embedding.

## Read in the data
First, let's read in the PBMCs data that we will use for the training. 

If you're using `colab` please make sure to upload your file(include `celltypes_PBMC.txt` and `counts_PBMC.csv`) first.
"""

expr_ctrl = pd.read_csv("/counts_PBMC.csv", sep = ",", index_col = 0).values
anno_ctrl = pd.read_csv("/celltypes_PBMC.txt", sep = "\t", header = None)
expr_ctrl = StandardScaler().fit_transform(expr_ctrl)

# Here we first reduce the dimension with PCA to speed up the training process
expr_ctrl = PCA(n_components = 100).fit_transform(expr_ctrl)
expr_ctrl = torch.FloatTensor(expr_ctrl)

"""## Create layers
The autoencoder that we use for biological dataset, e.g. single-cell dataset, is built upon multiple fully connected layers(fc layer). So our first step is to create our own fully connected layer from scratch.

**Note**: you may need to use [`nn.Parameter()`](https://pytorch.org/docs/stable/generated/torch.nn.parameter.Parameter.html?highlight=parameter).


"""

class fc_layer(nn.Module):
    def __init__(self, input_size, output_size, activation = nn.ReLU()):
        """\
          Parameters:
          ----------
          input_size: input dimension of fully connected layers
          output_size: output dimension of fully connected layers
          activation: activation function, default nn.ReLU()
        """
        super(fc_layer, self).__init__()

        # =========================================
        # weight matrix of the fully connected layer, please complete the code here
        # INSERT YOUR CODE BELOW
        # =========================================
        self.weight = nn.Parameter(torch.Tensor(input_size, output_size))
        # =========================================
        # bias of the fully connected layer
        # INSERT YOUR CODE BELOW
        # =========================================
        self.bias = nn.Parameter(torch.Tensor(1, output_size))
        # activation function of the fully connected layer
        self.activation = activation
        # initialize the parameters of weight and bias
        torch.nn.init.xavier_uniform_(self.weight)
        torch.nn.init.xavier_uniform_(self.bias)


    def forward(self, x):
        """\
          x: input data
        """
        # =========================================
        # implement forward pass here (before the activation function)
        # INSERT YOUR CODE BELOW
        # =========================================
        output = torch.mm(x, self.weight) + self.bias
        # call the activation function is there is one, you don't need to revise anything here
        if self.activation is not None:
          output = self.activation(output)
        return output

"""Next we create all the fc layers in the autoencoder using the `fc_layer` that we have already implemented above.

Layers for the autoencoder will be **number of genes** -> **128** -> **2** --> **128** -> **number of genes**. 

"""

# The first hidden layer has the input dimension the same as the number of genes in the single-cell dataset,
# the output dimension will be set as 128, relu will be used as the activation function. 
# =========================================
# Please complete the code for hidden_layer1 
# INSERT YOUR CODE BELOW
# =========================================
hidden_layer1 = fc_layer(100,128)

# The second hidden layer has the input dimension = 128, and the output dimension = 2. 
# No activation function will be used here.
# =========================================
# Please complete the code for hidden_layer2 
# INSERT YOUR CODE BELOW
# =========================================
hidden_layer2 = fc_layer(128,2, activation = None)

# The third hidden layer has the input dimension = 2, and the output dimension = 128. 
# Relu will be used as the activation function.
# =========================================
# Please complete the code for hidden_layer3
# INSERT YOUR CODE BELOW
# =========================================
hidden_layer3 = fc_layer(2,128)

# The output layer has the input dimension = 128, and the output dimension = number of genes. 
# No activation function will be used here.
# =========================================
# Please complete the code for output
# INSERT YOUR CODE BELOW
# =========================================
output = fc_layer(128,100)

"""## Assemble the autoencoder
After creating all the layers, the final step is to assemble the autoencoder with all those layers. `nn.Sequential` would be the most straightforward way fo the implementation. 

Please create your autoencoder using `nn.Sequential`, you might find the link below useful for your implementation:

* https://pytorch.org/docs/stable/generated/torch.nn.Sequential.html
* https://stackoverflow.com/questions/46141690/how-to-write-a-pytorch-sequential-model


"""

# =========================================
# Create the model using nn.Sequential
# INSERT YOUR CODE HERE
# =========================================
autoencoder = nn.Sequential(
 hidden_layer1,
 hidden_layer2,
 hidden_layer3,
 output,
)

# you can print the structure of the model out with print function
print(autoencoder)

"""## Define the optimizer
Gradient descent algorithm will be used to train the parameters of the autoencoder. After creating the autoencoder, we need to define the **optimizer**(can be considered as the gradient descent algorithm) for the model. 

`PyTorch` provides a number of different optimizers for us to choose from, such as `optim.SGD`, `optim.RMSprop`, `optim.Adagrad` and `optim.Adam`. Here we will use stochastic gradient descent(`optim.SGD`) for the optimization. 

"""

# Learning rate for the optimizer is a key hyperparameter to play around with when we need to train the model.
# Please try with different learning rate and find the one that is suitable for our task.
# =========================================
# INSERT YOUR CODE BELOW
# =========================================
learning_rate = 0.1

optimizer = optim.SGD(autoencoder.parameters(), lr=learning_rate)

"""## Construct loss function
Loss function is necessary when we want to train a DNN model. We usually use **Mean Square Error(MSE) loss** for autoencoder.

(**Note**: You don't need to implement any code here)
"""

loss_fcn = nn.MSELoss()

"""## Train the autoencoder
Here we have arrived at the final step, which is to train the model. We have already implemented this part for you. Make sure to carefully go through the code yourself, as it will help you gain a better understanding of how the model is trained.
"""

def train_model(model, optimizer, loss_fcn, n_epochs=10):

    batch_size=100
    losses = []

    # we'll train the network for 10 epochs
    step = 0
    for epoch in range(n_epochs):
        # randomize the order of the data each time through
        random_order = np.random.permutation(expr_ctrl.shape[0])
        data_randomized = expr_ctrl[random_order]

        # train the network on batches of size `batch_size`
        for data_batch in np.array_split(data_randomized, data_randomized.shape[0] // batch_size):
            step += 1

            # update the network weights to minimize the loss
            output = model(data_batch)

            # get loss
            loss = loss_fcn(output, data_batch)

            # print the loss every 100 epochs
            if step % 100 == 0:
                print("Step: {} Loss: {:.3f}".format(step, loss.item()))

            # backpropagate the loss
            loss.backward()

            # update parameters
            optimizer.step()

            # reset gradients
            optimizer.zero_grad()
            losses.append(loss)

    return losses

losses = train_model(autoencoder, optimizer = optimizer, loss_fcn = loss_fcn, n_epochs = 50)

"""## Plot the latent embedding
After training the autoencoder, we need to extract the latent embedding that it generates, which correspond to the output of the encoder. 
"""

# Get the output of the encoder(latent embedding) and visualize it with a scatter plot, complete the code below. Note that it is possible that the result might be different for different runs. 
with torch.no_grad():
    # note that ae_coordinates need to be transformed into numpy for the plot.
    # =========================================
    # INSERT YOUR CODE BELOW
    # =========================================
    ae_coordinates = autoencoder[:2](expr_ctrl).numpy()

def plot_latent(z, anno, save = None, figsize = (10,10), axis_label = "Latent", **kwargs):
    _kwargs = {
        "s": 10,
        "alpha": 0.9,
    }
    _kwargs.update(kwargs)

    fig = plt.figure(figsize = figsize)
    ax = fig.add_subplot()
    cluster_types = set([x for x in np.unique(anno)])
    colormap = plt.cm.get_cmap("tab20", len(cluster_types))

    for i, cluster_type in enumerate(cluster_types):
        index = np.where(anno == cluster_type)[0]
        ax.scatter(z[index,0], z[index,1], color = colormap(i), label = cluster_type, **_kwargs)
    
    ax.legend(loc='upper left', prop={'size': 15}, frameon = False, ncol = 1, bbox_to_anchor=(1.04, 1))
    
    ax.tick_params(axis = "both", which = "major", labelsize = 15)

    ax.set_xlabel(axis_label + " 1", fontsize = 19)
    ax.set_ylabel(axis_label + " 2", fontsize = 19)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)  
    
    if save:
        fig.savefig(save, bbox_inches = "tight")
    
    print(save)



plot_latent(ae_coordinates, anno_ctrl, axis_label = "Latent", save = "AE.pdf")

"""## Implement autoencoder using nn.Linear
`Pytorch` actually has wrapped up all the basic layer "Lego" for you, e.g. `nn.Linear` for fully connected layer, `nn.Conv2d` for convolutional layer, etc. 

Please implement the autoencoder we used above using `nn.Linear`, the backbone of the code is already provided as below:
"""

class autoencoder(nn.Module):
    def __init__(self, in_features):
        super(autoencoder,self).__init__()
        # complete the code for hidden_layer1, hidden_layer2, hidden_layer3, output, 
        # the dimensions of the layers are the same as above 
        # =========================================
        # INSERT YOUR CODE BELOW
        # =========================================
        self.hidden_layer1 = nn.Linear(100, 128)
        self.hidden_layer2 = nn.Linear(128, 2)
        self.hidden_layer3 = nn.Linear(2, 128)
        self.output = nn.Linear(128, 100)
    
    
    def forward(self, x):
        # complete the code for forward pass
        # =========================================
        # embed is the output of the encoder
        # INSERT YOUR CODE BELOW
        # =========================================
        embed = F.relu(self.hidden_layer1(x))
        embed = F.relu(self.hidden_layer2(embed))
        embed = F.relu(self.hidden_layer3(embed))
        # =========================================
        # output is the output of the decoder
        # INSERT YOUR CODE BELOW
        # =========================================
        output = torch.sigmoid(self.output(embed))


        return output

ae = autoencoder(in_features = expr_ctrl.shape[1])
optimizer = optim.SGD(ae.parameters(), lr=learning_rate)
losses = train_model(ae, optimizer = optimizer, loss_fcn = loss_fcn, n_epochs = 100)

# Get the output of the encoder(latent embedding) and visualize it with a scatter plot, complete the code below. Note that plot that you get might not be the same as the first plot, and it's also possible that the result is different for different runs.
with torch.no_grad():
    # note that ae_coordinates need to be transformed into numpy for the plot.
    # =========================================
    # INSERT YOUR CODE BELOW
    # =========================================
    ae_coordinates = ae.hidden_layer2(F.relu(ae.hidden_layer1(expr_ctrl)))

plot_latent(ae_coordinates, anno_ctrl, axis_label = "Latent", save = "AE2.pdf")